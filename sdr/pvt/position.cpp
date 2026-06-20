#include "position.h"
#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include "logging.h"

namespace
{
// Unit line-of-sight vectors (receiver -> each satellite) and their ranges.
void los_vector( const Eigen::MatrixXd& xusr, const Eigen::MatrixXd& xsat, Eigen::MatrixXd& u, Eigen::VectorXd& rng )
{
    u   = xsat.rowwise() - xusr.row( 0 );
    rng = u.rowwise().norm();
    for( int i = 0; i < u.rows(); ++i )
    {
        u.row( i ) /= rng( i );
    }
}

// Geometry matrix for a position+clock solve: rows [-u, 1].
Eigen::MatrixXd jacobian_residuals( const Eigen::VectorXd& x, const Eigen::MatrixXd& xsat )
{
    Eigen::MatrixXd u;
    Eigen::VectorXd rng;
    los_vector( x.head( 3 ).transpose(), xsat, u, rng );

    Eigen::MatrixXd J( u.rows(), 4 );
    J.leftCols( 3 ) = -u;
    J.col( 3 )      = Eigen::VectorXd::Ones( u.rows() );
    return J;
}

Eigen::VectorXd pr_residuals( const Eigen::VectorXd& x, const Eigen::MatrixXd& xsat, const Eigen::VectorXd& pr )
{
    Eigen::MatrixXd u;
    Eigen::VectorXd rng;
    los_vector( x.head( 3 ).transpose(), xsat, u, rng );
    return rng - ( pr.array() - x( 3 ) ).matrix();
}

Eigen::VectorXd prr_residuals(
    const Eigen::VectorXd& v,
    const Eigen::MatrixXd& vsat,
    const Eigen::VectorXd& prr,
    const Eigen::VectorXd& x,
    const Eigen::MatrixXd& xsat
)
{
    Eigen::MatrixXd u;
    Eigen::VectorXd rng;
    los_vector( x.head( 3 ).transpose(), xsat, u, rng );

    Eigen::VectorXd rate( u.rows() );
    for( int i = 0; i < u.rows(); ++i )
    {
        rate( i ) = ( vsat.row( i ).head( 3 ) - v.head( 3 ).transpose() ).dot( u.row( i ) );
    }
    return rate - ( prr.array() - v( 3 ) ).matrix();
}

Eigen::VectorXd least_squares(
    const Eigen::VectorXd&                                          v0,
    const std::function<Eigen::VectorXd( const Eigen::VectorXd& )>& residuals_func,
    const std::function<Eigen::MatrixXd( const Eigen::VectorXd& )>& jacobian_func,
    int                                                             max_iters,
    double                                                          tol
)
{
    Eigen::VectorXd v = v0;
    for( int iter = 0; iter < max_iters; ++iter )
    {
        Eigen::VectorXd r     = residuals_func( v );
        Eigen::MatrixXd J     = jacobian_func( v );
        Eigen::VectorXd delta = -( J.transpose() * J ).ldlt().solve( J.transpose() * r );
        v += delta;
        if( delta.norm() < tol )
        {
            break;
        }
    }
    return v;
}

template <typename T, typename MemberType>
Eigen::Matrix<MemberType, Eigen::Dynamic, 1> extract_vector( const std::vector<T>& data, MemberType T::* member_ptr )
{
    static_assert( std::is_arithmetic<MemberType>::value, "Member must be scalar" );
    Eigen::Matrix<MemberType, Eigen::Dynamic, 1> vec( data.size() );
    for( size_t i = 0; i < data.size(); ++i )
    {
        vec( i ) = data[i].*member_ptr;
    }
    return vec;
}

// Pack satellite positions / velocities into an n x 3 matrix.
Eigen::MatrixXd sat_matrix(
    const std::vector<Satellite_measurement>& m,
    double Satellite_measurement::* px,
    double Satellite_measurement::* py,
    double Satellite_measurement::* pz
)
{
    Eigen::MatrixXd out( m.size(), 3 );
    for( size_t i = 0; i < m.size(); ++i )
    {
        out( i, 0 ) = m[i].*px;
        out( i, 1 ) = m[i].*py;
        out( i, 2 ) = m[i].*pz;
    }
    return out;
}

// State index layout: [0..2] position, [3..5] velocity, [6] clock bias, [7] clock drift,
// [Position_solver::ISB_BASE + constellation] inter-system bias (one slot per constellation).
constexpr int PX = 0, VX = 3, CB = 6, CD = 7;

// Mixed-constellation seed LS. State x = [px py pz cb isb_0 isb_1 ...]; isb_col[i] is measurement
// i's inter-system-bias column (>=4), or -1 for the reference constellation (no isb term). This is
// what lets the filter seed off a MIX of constellations (position + clock + one bias per non-ref).
Eigen::VectorXd pr_residuals_isb(
    const Eigen::VectorXd& x, const Eigen::MatrixXd& xsat, const Eigen::VectorXd& pr, const std::vector<int>& isb_col
)
{
    Eigen::MatrixXd u;
    Eigen::VectorXd rng;
    los_vector( x.head( 3 ).transpose(), xsat, u, rng );
    Eigen::VectorXd res( xsat.rows() );
    for( int i = 0; i < xsat.rows(); ++i )
    {
        const double clk = x( 3 ) + ( isb_col[i] >= 0 ? x( isb_col[i] ) : 0.0 );
        res( i )         = rng( i ) - ( pr( i ) - clk );
    }
    return res;
}

Eigen::MatrixXd
jacobian_isb( const Eigen::VectorXd& x, const Eigen::MatrixXd& xsat, const std::vector<int>& isb_col, int nstate )
{
    Eigen::MatrixXd u;
    Eigen::VectorXd rng;
    los_vector( x.head( 3 ).transpose(), xsat, u, rng );
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero( xsat.rows(), nstate );
    J.leftCols( 3 )   = -u;
    J.col( 3 ).setOnes();
    for( int i = 0; i < xsat.rows(); ++i )
    {
        if( isb_col[i] >= 0 )
        {
            J( i, isb_col[i] ) = 1.0;
        }
    }
    return J;
}

// Robust velocity + clock-drift (RAIM). The velocity/clock-drift split from the pseudorange rates
// is fragile: a CLUSTER of correlated outliers (e.g. several marginally-locked Galileo SVs that all
// agree with each other) defeats plain largest-residual rejection - it locks onto the bad cluster's
// consensus and throws out the good SVs instead. So we ANCHOR the velocity on the reference
// constellation (the trusted majority that seeded the fix, with good geometry), then accept every
// other SV's prr only if it agrees with that anchored velocity (residual <= max_residual). Mirrors
// the reference-constellation logic used for position/clock. Returns the per-measurement inlier mask
// (only the prr / velocity is screened; the pr / position rows are untouched) and fills vel_out.
std::vector<bool> robust_velocity(
    const Eigen::MatrixXd&   sat_vel,
    const Eigen::MatrixXd&   sat_pos,
    const Eigen::VectorXd&   prr,
    const Eigen::VectorXd&   pos,
    const std::vector<bool>& is_reference,
    double                   max_residual,
    Eigen::VectorXd&         vel_out
)
{
    const int n = static_cast<int>( prr.size() );

    // Solve [vx vy vz cd] from a subset of SV indices.
    auto solve = [&]( const std::vector<int>& idx ) -> Eigen::VectorXd
    {
        const int       m = static_cast<int>( idx.size() );
        Eigen::MatrixXd sv( m, 3 ), sp( m, 3 );
        Eigen::VectorXd pr_in( m );
        for( int k = 0; k < m; ++k )
        {
            sv.row( k ) = sat_vel.row( idx[k] );
            sp.row( k ) = sat_pos.row( idx[k] );
            pr_in( k )  = prr( idx[k] );
        }
        return least_squares(
            Eigen::VectorXd::Zero( 4 ),
            [&]( const Eigen::VectorXd& vv ) { return prr_residuals( vv, sv, pr_in, pos, sp ); },
            [&]( const Eigen::VectorXd& ) { return jacobian_residuals( pos, sp ); },
            50,
            1e-6
        );
    };

    std::vector<int> ref_idx, all_idx;
    for( int i = 0; i < n; ++i )
    {
        all_idx.push_back( i );
        if( is_reference[i] )
        {
            ref_idx.push_back( i );
        }
    }

    // Anchor on the reference constellation if it can solve on its own; else fall back to all SVs.
    vel_out = solve( static_cast<int>( ref_idx.size() ) >= 4 ? ref_idx : all_idx );
    if( vel_out.hasNaN() )
    {
        return std::vector<bool>( n, true );
    }

    // Inliers = SVs whose prr agrees with the anchored velocity (|prr_residual| == |innovation|).
    Eigen::VectorXd   r = prr_residuals( vel_out, sat_vel, prr, pos, sat_pos );
    std::vector<bool> inlier( n );
    std::vector<int>  inlier_idx;
    for( int i = 0; i < n; ++i )
    {
        inlier[i] = std::abs( r( i ) ) <= max_residual;
        if( inlier[i] )
        {
            inlier_idx.push_back( i );
        }
    }

    // Refine over all inliers for accuracy (the anchor stays the trust source).
    if( static_cast<int>( inlier_idx.size() ) >= 4 )
    {
        Eigen::VectorXd refined = solve( inlier_idx );
        if( !refined.hasNaN() )
        {
            vel_out = refined;
        }
    }
    return inlier;
}

// Highest-priority (enum order) constellation present (>=1 measurement). Used as the seed gauge:
// its inter-system bias is the zero, every other constellation's bias is measured relative to it.
bool pick_reference( const std::vector<Satellite_measurement>& meas, int min_sats, Constellation& out )
{
    int count[NUM_CONSTELLATIONS] = { 0 };
    for( const auto& s : meas )
    {
        const int idx = static_cast<int>( s.constellation );
        if( idx >= 0 && idx < NUM_CONSTELLATIONS )
        {
            ++count[idx];
        }
    }
    for( int c = 0; c < NUM_CONSTELLATIONS; ++c )
    {
        if( count[c] >= min_sats )
        {
            out = static_cast<Constellation>( c );
            return true;
        }
    }
    return false;
}

} // namespace

bool Position_solver::initialise( const std::vector<Satellite_measurement>& measurements )
{
    if( measurements.size() < 4 )
    {
        return false;
    }

    // Reference constellation = highest-priority one present; its inter-system bias is the gauge
    // zero (folded into the clock bias). Each OTHER constellation present gets its own isb column,
    // so the seed spans a MIX of constellations (position + clock + one bias per non-reference).
    if( !pick_reference( measurements, 1, reference_ ) )
    {
        return false;
    }

    int ls_col[NUM_CONSTELLATIONS];
    for( int c = 0; c < NUM_CONSTELLATIONS; ++c )
    {
        ls_col[c] = -1;
    }

    int nstate = 4; // [px py pz cb] + one column per non-reference constellation present
    {
        bool present[NUM_CONSTELLATIONS] = { false };
        for( const auto& m : measurements )
        {
            present[static_cast<int>( m.constellation )] = true;
        }
        for( int c = 0; c < NUM_CONSTELLATIONS; ++c )
        {
            if( present[c] && static_cast<Constellation>( c ) != reference_ )
            {
                ls_col[c] = nstate++;
            }
        }
    }

    if( static_cast<int>( measurements.size() ) < nstate )
    {
        return false; // need >= 4 + (#non-reference constellations) sats to solve
    }

    // Seed pseudorange RAIM. A gross outlier - e.g. a cross-correlation / false-preamble track whose
    // garbage range would drag the seed hundreds of km off (seen on the Skydel sim, where one false-locked
    // SV pushed the seed to |pos| ~8380 km) - is NOT caught by the PDOP gate (geometry is fine, the
    // measurement is not) and then poisons the whole filter. We can't use the reference-anchored trick the
    // velocity RAIM uses, because the bad SV may be IN the reference constellation (the Galileo-only case).
    // We also can't just drop the largest post-fit residual: a high-LEVERAGE outlier smears its error onto
    // the inliers, so the biggest residual can land on a GOOD SV. So use LEAVE-ONE-OUT: while the fit's max
    // residual is grossly large, drop the SV whose REMOVAL most reduces that max - which cleanly identifies
    // a single outlier regardless of leverage - and re-fit. A drop is allowed only while it leaves >= nstate
    // measurements AND keeps every present constellation with a satellite (nstate / ls_col are fixed from the
    // full set, so the reference gauge + each inter-system-bias column stay observable). A coordinated
    // CLUSTER of outliers could still fool this, but the PDOP gate + the EKF innovation gate guard that.
    auto fit = [&]( const std::vector<Satellite_measurement>& w, Eigen::VectorXd* pos_out,
                    Eigen::MatrixXd* sat_pos_out, std::vector<int>* col_out ) -> double
    {
        Eigen::MatrixXd sp = sat_matrix(
            w,
            &Satellite_measurement::satellite_pos_x,
            &Satellite_measurement::satellite_pos_y,
            &Satellite_measurement::satellite_pos_z
        );
        Eigen::VectorXd  pr = extract_vector( w, &Satellite_measurement::pseudorange_m );
        std::vector<int> c( w.size() );
        for( size_t i = 0; i < w.size(); ++i )
        {
            c[i] = ls_col[static_cast<int>( w[i].constellation )];
        }
        Eigen::VectorXd p = least_squares(
            Eigen::VectorXd::Zero( nstate ),
            [&]( const Eigen::VectorXd& x ) { return pr_residuals_isb( x, sp, pr, c ); },
            [&]( const Eigen::VectorXd& x ) { return jacobian_isb( x, sp, c, nstate ); },
            50,
            1e-6
        );
        if( p.hasNaN() )
        {
            return -1.0;
        }
        const Eigen::VectorXd r  = pr_residuals_isb( p, sp, pr, c );
        double                mx = 0.0;
        for( int i = 0; i < r.size(); ++i )
        {
            mx = std::max( mx, std::abs( r( i ) ) );
        }
        if( pos_out )
            *pos_out = p;
        if( sat_pos_out )
            *sat_pos_out = sp;
        if( col_out )
            *col_out = c;
        return mx;
    };

    std::vector<Satellite_measurement> work = measurements;
    Eigen::VectorXd                    pos;
    Eigen::MatrixXd                    sat_pos;
    std::vector<int>                   col;
    for( ;; )
    {
        const double mx = fit( work, &pos, &sat_pos, &col );
        if( mx < 0.0 )
        {
            return false; // diverged
        }
        if( mx <= SEED_PR_RAIM_RESIDUAL_M || static_cast<int>( work.size() ) <= nstate )
        {
            break; // good fit, or can't drop any more
        }

        int cc[NUM_CONSTELLATIONS] = { 0 };
        for( const auto& m : work )
        {
            cc[static_cast<int>( m.constellation )]++;
        }

        int    best_i   = -1;
        double best_max = mx; // a removal must strictly improve the max residual
        for( int i = 0; i < static_cast<int>( work.size() ); ++i )
        {
            if( cc[static_cast<int>( work[i].constellation )] < 2 )
            {
                continue; // keep every present constellation populated
            }
            std::vector<Satellite_measurement> trial = work;
            trial.erase( trial.begin() + i );
            const double m2 = fit( trial, nullptr, nullptr, nullptr );
            if( m2 >= 0.0 && m2 < best_max )
            {
                best_max = m2;
                best_i   = i;
            }
        }
        if( best_i < 0 )
        {
            break; // no single removal helps - leave it to the PDOP / innovation gates
        }
        work.erase( work.begin() + best_i );
    }

    // Geometry gate: reject a poorly-determined (high-PDOP) seed - e.g. a sparse single-
    // constellation set with bad vertical geometry - rather than locking the filter onto a bad
    // fix. This REPLACES privileging a constellation: any mix with good geometry may seed, and we
    // simply wait for more/better sats otherwise (no dependence on which constellation, or order).
    Eigen::MatrixXd J    = jacobian_isb( pos, sat_pos, col, nstate );
    Eigen::MatrixXd cov  = ( J.transpose() * J ).inverse();
    double          pdop = std::sqrt( cov( 0, 0 ) + cov( 1, 1 ) + cov( 2, 2 ) );

    if( !std::isfinite( pdop ) || pdop > MAX_SEED_PDOP )
    {
        return false;
    }

    // Velocity + clock-drift via robust (RAIM) LS - the clock drift is shared (one oscillator), so
    // there is no per-constellation drift term. RAIM anchors on the reference constellation and drops
    // SVs whose Doppler disagrees, so a biased-Doppler SV can't corrupt the velocity seed (which the
    // EKF would otherwise lock onto). Built from the SAME inlier set the position RAIM accepted.
    Eigen::MatrixXd sat_vel = sat_matrix(
        work,
        &Satellite_measurement::satellite_vel_x,
        &Satellite_measurement::satellite_vel_y,
        &Satellite_measurement::satellite_vel_z
    );
    Eigen::VectorXd   prr = extract_vector( work, &Satellite_measurement::pseudorange_rate_m_s );
    std::vector<bool> is_ref( work.size() );
    for( size_t i = 0; i < work.size(); ++i )
    {
        is_ref[i] = ( work[i].constellation == reference_ );
    }

    Eigen::VectorXd vel;
    robust_velocity( sat_vel, sat_pos, prr, pos, is_ref, VELOCITY_RAIM_RESIDUAL_M_S, vel );

    if( vel.hasNaN() )
    {
        return false;
    }

    x_.setZero();
    x_.segment<3>( PX ) = pos.head<3>();
    x_.segment<3>( VX ) = vel.head<3>();
    x_( CB )            = pos( 3 );
    x_( CD )            = vel( 3 );

    for( int c = 0; c < NUM_CONSTELLATIONS; ++c )
    {
        if( ls_col[c] >= 0 )
        {
            x_( ISB_BASE + c ) = pos( ls_col[c] );
        }
    }

    // Initial covariance: generous, so the first few updates pull the state in quickly.
    P_.setZero();
    P_( PX + 0, PX + 0 ) = P_( PX + 1, PX + 1 ) = P_( PX + 2, PX + 2 ) = 900.0; // position    (30 m)^2
    P_( VX + 0, VX + 0 ) = P_( VX + 1, VX + 1 ) = P_( VX + 2, VX + 2 ) = 100.0; // velocity    (10 m/s)^2
    P_( CB, CB )                                                       = 1.0e6; // clock bias   (1 km)^2
    P_( CD, CD )                                                       = 1.0e4; // clock drift  (100 m/s)^2

    for( int c = 0; c < NUM_CONSTELLATIONS; ++c )
    {
        if( static_cast<Constellation>( c ) == reference_ )
        {
            continue; // reference bias pinned at 0 (P=0, no process noise)
        }
        // A bias estimated in THIS seed starts moderately known; one for a constellation not yet
        // seen starts wholly unknown (huge), so the gate accepts its first (ms-level) update later.
        P_( ISB_BASE + c, ISB_BASE + c ) = ( ls_col[c] >= 0 ) ? 1.0e6 : 1.0e12;
    }

    // Fresh seed: clear the latched present-constellation set, then record this seed's mix.
    std::fill( std::begin( isb_present_ ), std::end( isb_present_ ), false );
    latch_present( measurements );

    initialised_ = true;
    return true;
}

void Position_solver::predict( double dt )
{
    // State transition F: position integrates velocity (p += v*dt) and clock bias integrates
    // clock drift (cb += cd*dt). Velocity and drift themselves hold.
    Eigen::Matrix<double, N, N> F = Eigen::Matrix<double, N, N>::Identity();
    for( int i = 0; i < 3; ++i )
    {
        F( PX + i, VX + i ) = dt;
    }
    F( CB, CD ) = dt;

    // Process noise Q. Spatial part: per-axis white-noise-acceleration over each [p, v] pair.
    // Clock part: the standard [bias, drift] block - bias phase noise (CLOCK_BIAS_PSD) plus the
    // drift's contribution integrated into the bias (CLOCK_DRIFT_PSD terms).
    Eigen::Matrix<double, N, N> Q = Eigen::Matrix<double, N, N>::Zero();

    const double dt2 = dt * dt, dt3 = dt2 * dt;
    const double q_pp = ACCEL_PSD * dt3 / 3.0;
    const double q_pv = ACCEL_PSD * dt2 / 2.0;
    const double q_vv = ACCEL_PSD * dt;

    for( int i = 0; i < 3; ++i )
    {
        Q( PX + i, PX + i ) = q_pp;
        Q( PX + i, VX + i ) = q_pv;
        Q( VX + i, PX + i ) = q_pv;
        Q( VX + i, VX + i ) = q_vv;
    }

    Q( CB, CB ) = CLOCK_BIAS_PSD * dt + CLOCK_DRIFT_PSD * dt3 / 3.0;
    Q( CB, CD ) = CLOCK_DRIFT_PSD * dt2 / 2.0;
    Q( CD, CB ) = CLOCK_DRIFT_PSD * dt2 / 2.0;
    Q( CD, CD ) = CLOCK_DRIFT_PSD * dt;

    // Inter-system biases: slow random walks (F holds them via identity). Only non-reference
    // slots get process noise; the reference slot stays pinned at 0.
    for( int c = 0; c < NUM_CONSTELLATIONS; ++c )
    {
        if( static_cast<Constellation>( c ) != reference_ )
        {
            Q( ISB_BASE + c, ISB_BASE + c ) = ISB_PSD * dt;
        }
    }

    x_ = F * x_;
    P_ = F * P_ * F.transpose() + Q;
}

void Position_solver::update( const std::vector<Satellite_measurement>& measurements )
{
    const int n = static_cast<int>( measurements.size() );

    const Eigen::Vector3d p  = x_.segment<3>( PX );
    const Eigen::Vector3d v  = x_.segment<3>( VX );
    const double          cb = x_( CB );
    const double          cd = x_( CD );

    // Velocity RAIM: flag which pseudorange-rates are mutually consistent (a consensus velocity solve),
    // INDEPENDENT of the current EKF velocity. A biased-Doppler SV is then excluded from the prr update
    // every tick, so one bad SV can't pull the filter onto - or hold it at - a wrong velocity (the pr /
    // position rows are untouched: that SV's code lock is usually fine).
    Eigen::VectorXd prr_v     = extract_vector( measurements, &Satellite_measurement::pseudorange_rate_m_s );
    Eigen::MatrixXd sat_pos_m = sat_matrix(
        measurements,
        &Satellite_measurement::satellite_pos_x,
        &Satellite_measurement::satellite_pos_y,
        &Satellite_measurement::satellite_pos_z
    );
    Eigen::MatrixXd sat_vel_m = sat_matrix(
        measurements,
        &Satellite_measurement::satellite_vel_x,
        &Satellite_measurement::satellite_vel_y,
        &Satellite_measurement::satellite_vel_z
    );
    Eigen::VectorXd   pos3 = x_.segment<3>( PX );
    std::vector<bool> is_ref( n );
    for( int i = 0; i < n; ++i )
    {
        is_ref[i] = ( measurements[i].constellation == reference_ );
    }
    Eigen::VectorXd   raim_vel;
    std::vector<bool> prr_inlier =
        robust_velocity( sat_vel_m, sat_pos_m, prr_v, pos3, is_ref, VELOCITY_RAIM_RESIDUAL_M_S, raim_vel );

    // Build each candidate measurement row, then keep only those whose innovation is plausible
    // (normalised innovation^2 = y^2 / (H P H^T + R) within the chi-square gate). This rejects
    // outliers - stale/garbage observables as channels drain at end-of-file, gross multipath -
    // so a few bad measurements cannot drag the filter off. The same row is also weighted by
    // satellite elevation through R.
    std::vector<Eigen::Matrix<double, 1, N>> rows;
    std::vector<double>                      ys;
    std::vector<double>                      Rs;
    rows.reserve( 2 * n );
    ys.reserve( 2 * n );
    Rs.reserve( 2 * n );

    auto consider = [&]( const Eigen::Matrix<double, 1, N>& h, double innovation, double variance )
    {
        const double s = h * P_ * h.transpose() + variance;
        if( innovation * innovation <= INNOV_GATE2 * s )
        {
            rows.push_back( h );
            ys.push_back( innovation );
            Rs.push_back( variance );
        }
    };

    for( int i = 0; i < n; ++i )
    {
        const Satellite_measurement& sm = measurements[i];
        Eigen::Vector3d              sat_p( sm.satellite_pos_x, sm.satellite_pos_y, sm.satellite_pos_z );
        Eigen::Vector3d              sat_v( sm.satellite_vel_x, sm.satellite_vel_y, sm.satellite_vel_z );

        Eigen::Vector3d los   = sat_p - p;
        double          range = los.norm();
        Eigen::Vector3d u     = los / range; // unit receiver->satellite

        // Elevation weighting: variance scales as 1/sin^2(el). Unknown elevation (0, before a
        // fix has fed the observation engine) falls back to the base variance.
        double sin_el = 1.0;
        if( sm.elevation_rad > 0.0 )
        {
            sin_el = std::sin( std::max( sm.elevation_rad, MIN_ELEVATION_RAD ) );
        }
        const double w = 1.0 / ( sin_el * sin_el );

        // Pseudorange row: pr = range + clock_bias (+ inter-system bias for non-reference
        // constellations). d/dp = -u, d/dcb = 1, d/disb[c] = 1 (non-reference only).
        Eigen::Matrix<double, 1, N> h_pr = Eigen::Matrix<double, 1, N>::Zero();
        h_pr.segment<3>( PX )            = -u.transpose();
        h_pr( CB )                       = 1.0;
        double pr_pred                   = range + cb;
        if( sm.constellation != reference_ )
        {
            const int isb_i = ISB_BASE + static_cast<int>( sm.constellation );
            h_pr( isb_i )   = 1.0;
            pr_pred += x_( isb_i );
        }
        consider( h_pr, sm.pseudorange_m - pr_pred, PR_STD_M * PR_STD_M * w );

        // Pseudorange-rate row: prr = (vsat - v).u + clock_drift.  d/dv = -u, d/dcd = 1.
        // (The small d/dp term from u rotating is negligible and dropped, as is standard.)
        // Skipped for an SV the velocity RAIM flagged as a biased-Doppler outlier.
        if( prr_inlier[i] )
        {
            Eigen::Matrix<double, 1, N> h_prr = Eigen::Matrix<double, 1, N>::Zero();
            h_prr.segment<3>( VX )            = -u.transpose();
            h_prr( CD )                       = 1.0;
            consider( h_prr, sm.pseudorange_rate_m_s - ( ( sat_v - v ).dot( u ) + cd ), PRR_STD_M_S * PRR_STD_M_S * w );
        }
    }

    const int m = static_cast<int>( rows.size() );
    if( m == 0 )
    {
        return; // nothing passed the gate - coast on the prediction
    }

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero( m, N );
    Eigen::VectorXd y( m );
    Eigen::VectorXd R( m );
    for( int i = 0; i < m; ++i )
    {
        H.row( i ) = rows[i];
        y( i )     = ys[i];
        R( i )     = Rs[i];
    }

    // Standard EKF update. S = H P H^T + R; K = P H^T S^-1.
    Eigen::MatrixXd Ht = H.transpose();
    Eigen::MatrixXd S  = H * P_ * Ht;
    S.diagonal() += R;
    Eigen::MatrixXd K = P_ * Ht * S.inverse();

    x_ += K * y;
    Eigen::Matrix<double, N, N> I = Eigen::Matrix<double, N, N>::Identity();
    P_                            = ( I - K * H ) * P_;
}

Position_solution Position_solver::as_solution() const
{
    Position_solution out;
    out.valid           = true;
    out.ecef_x_m        = x_( PX + 0 );
    out.ecef_y_m        = x_( PX + 1 );
    out.ecef_z_m        = x_( PX + 2 );
    out.clock_bias_m    = x_( CB );
    out.ecef_x_m_s      = x_( VX + 0 );
    out.ecef_y_m_s      = x_( VX + 1 );
    out.ecef_z_m_s      = x_( VX + 2 );
    out.clock_drift_m_s = x_( CD );

    // Reference constellation + the per-constellation inter-system biases (shown only for the
    // non-reference constellations that have actually contributed).
    out.reference = reference_;
    for( int c = 0; c < NUM_CONSTELLATIONS; ++c )
    {
        out.isb_m[c]       = x_( ISB_BASE + c );
        out.isb_present[c] = isb_present_[c];
    }
    return out;
}

void Position_solver::latch_present( const std::vector<Satellite_measurement>& measurements )
{
    for( const Satellite_measurement& m : measurements )
    {
        const int c = static_cast<int>( m.constellation );
        if( c >= 0 && c < NUM_CONSTELLATIONS && m.constellation != reference_ )
        {
            isb_present_[c] = true;
        }
    }
}

std::optional<Position_solution>
Position_solver::compute_solution( const std::vector<Satellite_measurement>& measurements, double rx_time_s )
{
    if( !initialised_ )
    {
        if( measurements.size() < 4 || !initialise( measurements ) )
        {
            return std::nullopt;
        }
        last_time_s_ = rx_time_s;
        return as_solution();
    }

    // Predict forward to this epoch. Guard dt against a stalled / non-monotonic clock.
    double dt = rx_time_s - last_time_s_;
    if( dt <= 0.0 || dt > 10.0 )
    {
        dt = 0.0;
    }
    predict( dt );
    last_time_s_ = rx_time_s;

    // Fuse the measurements (a full update needs >=4 sats; otherwise coast on the prediction).
    if( measurements.size() >= 4 )
    {
        update( measurements );
        latch_present( measurements );
    }

    if( !x_.allFinite() )
    {
        logging::log( logging::Level::Warning, "EKF state went non-finite; reinitialising on the next epoch" );
        initialised_ = false;
        return std::nullopt;
    }
    return as_solution();
}

#ifdef ENABLE_UNIT_TESTS
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

namespace
{
const double DEG = M_PI / 180.0;

Ecef geodetic_to_ecef( double lat, double lon, double alt )
{
    constexpr double A = 6378137.0, F = 1.0 / 298.257223563, E2 = F * ( 2.0 - F );
    const double     s = std::sin( lat );
    const double     N = A / std::sqrt( 1.0 - E2 * s * s );
    return { ( N + alt ) * std::cos( lat ) * std::cos( lon ), ( N + alt ) * std::cos( lat ) * std::sin( lon ),
             ( N * ( 1.0 - E2 ) + alt ) * s };
}

// A satellite at azimuth/elevation (rad) and slant range from the user, placed in ECEF via the
// local ENU frame at the user's geodetic location - guarantees a known elevation + good geometry.
Ecef sat_from_look( const Ecef& user, double lat, double lon, double az, double el, double range )
{
    const double e = std::cos( el ) * std::sin( az ), n = std::cos( el ) * std::cos( az ), u = std::sin( el );
    const double sl = std::sin( lat ), cl = std::cos( lat ), so = std::sin( lon ), co = std::cos( lon );
    const double dx = -so * e - sl * co * n + cl * co * u;
    const double dy = co * e - sl * so * n + cl * so * u;
    const double dz = cl * n + sl * u;
    return { user.x + range * dx, user.y + range * dy, user.z + range * dz };
}
} // namespace

TEST_CASE( "ekf_converges_to_known_truth", "[pvt][position][ekf]" )
{
    // Truth: a static receiver with a known clock bias + drift, and GPS satellites spread across the
    // sky. Build NOISE-FREE pseudoranges pr = |sat - user| + clock_bias (the sat-clock/iono/Sagnac
    // corrections live upstream in observation.cpp, not the solver), and rates prr = clock_drift
    // (static rx, static sats). The filter must recover the truth exactly (to numerical precision).
    const double lat = 45.0 * DEG, lon = 10.0 * DEG, alt = 120.0;
    const Ecef   user = geodetic_to_ecef( lat, lon, alt );
    const double cb_truth = 12345.0;   // receiver clock bias (m)
    const double cd_truth = -1600.0;   // receiver clock drift (m/s)

    struct
    {
        double az_deg, el_deg;
    } look[] = { { 30, 72 }, { 95, 18 }, { 150, 45 }, { 210, 28 }, { 275, 61 }, { 340, 35 }, { 60, 12 } };

    // Static satellites: fixed ECEF positions + geometric ranges to the (static) user.
    std::vector<Ecef>   sat_pos;
    std::vector<double> range;
    for( const auto& lk : look )
    {
        const Ecef s = sat_from_look( user, lat, lon, lk.az_deg * DEG, lk.el_deg * DEG, 22.0e6 );
        sat_pos.push_back( s );
        range.push_back( std::sqrt(
            ( s.x - user.x ) * ( s.x - user.x ) + ( s.y - user.y ) * ( s.y - user.y ) + ( s.z - user.z ) * ( s.z - user.z )
        ) );
    }

    Position_solver   solver;
    Position_solution sol {};
    const int         epochs = 40;
    double            t      = 0.0;
    for( int epoch = 0; epoch < epochs; ++epoch )
    {
        t = static_cast<double>( epoch );
        // Time-consistent measurements: the clock bias evolves as cb(t) = cb0 + cd*t, so the
        // pseudorange does too (the coupled EKF clock model expects exactly this). The rate is cd.
        const double cb_t = cb_truth + cd_truth * t;
        std::vector<Satellite_measurement> meas;
        for( size_t i = 0; i < sat_pos.size(); ++i )
        {
            Satellite_measurement m {};
            m.satellite_id         = static_cast<Satellite_id>( i + 1 );
            m.constellation        = Constellation::Gps;
            m.satellite_pos_x      = sat_pos[i].x;
            m.satellite_pos_y      = sat_pos[i].y;
            m.satellite_pos_z      = sat_pos[i].z;
            m.pseudorange_m        = range[i] + cb_t;
            m.pseudorange_rate_m_s = cd_truth; // static sat + static rx -> rate is pure clock drift
            m.elevation_rad        = look[i].el_deg * DEG;
            meas.push_back( m );
        }
        auto out = solver.compute_solution( meas, t );
        REQUIRE( out.has_value() ); // seeds on the first epoch (good geometry, 7 sats)
        sol = *out;
    }

    REQUIRE( sol.valid );
    REQUIRE( sol.reference == Constellation::Gps );
    REQUIRE_FALSE( sol.isb_present[static_cast<int>( Constellation::Gps )] ); // single-constellation: no ISB
    REQUIRE( sol.ecef_x_m == Catch::Approx( user.x ).margin( 1e-3 ) );
    REQUIRE( sol.ecef_y_m == Catch::Approx( user.y ).margin( 1e-3 ) );
    REQUIRE( sol.ecef_z_m == Catch::Approx( user.z ).margin( 1e-3 ) );
    REQUIRE( sol.clock_bias_m == Catch::Approx( cb_truth + cd_truth * t ).margin( 1e-2 ) );
    // Velocity/drift states converge over the run to the static truth.
    REQUIRE( sol.ecef_x_m_s == Catch::Approx( 0.0 ).margin( 1e-2 ) );
    REQUIRE( sol.ecef_y_m_s == Catch::Approx( 0.0 ).margin( 1e-2 ) );
    REQUIRE( sol.ecef_z_m_s == Catch::Approx( 0.0 ).margin( 1e-2 ) );
    REQUIRE( sol.clock_drift_m_s == Catch::Approx( cd_truth ).margin( 1e-2 ) );
}

TEST_CASE( "seed_position_raim_rejects_gross_outlier", "[pvt][position][raim]" )
{
    // Six clean GPS SVs (good geometry) consistent with a known rx position + clock, PLUS one false-lock
    // SV whose pseudorange is grossly wrong (off by 5000 km). Without seed RAIM the single outlier drags
    // the least-squares seed hundreds of km off (the Skydel Galileo-alone failure); with the largest-
    // residual RAIM in initialise() the outlier is dropped and the seed recovers truth from the 6 good SVs.
    const double lat = 60.18 * DEG, lon = 24.83 * DEG, alt = 47.0; // ~Helsinki (the Skydel truth)
    const Ecef   user     = geodetic_to_ecef( lat, lon, alt );
    const double cb_truth = 8000.0;

    struct
    {
        double az_deg, el_deg;
    } look[] = { { 30, 72 }, { 95, 18 }, { 150, 45 }, { 210, 28 }, { 275, 61 }, { 340, 35 }, { 60, 15 } };

    std::vector<Satellite_measurement> meas;
    for( int i = 0; i < 7; ++i )
    {
        const Ecef   s   = sat_from_look( user, lat, lon, look[i].az_deg * DEG, look[i].el_deg * DEG, 22.0e6 );
        const double rng = std::sqrt(
            ( s.x - user.x ) * ( s.x - user.x ) + ( s.y - user.y ) * ( s.y - user.y )
            + ( s.z - user.z ) * ( s.z - user.z )
        );
        Satellite_measurement m {};
        m.satellite_id         = static_cast<Satellite_id>( i + 1 );
        m.constellation        = Constellation::Gps;
        m.satellite_pos_x      = s.x;
        m.satellite_pos_y      = s.y;
        m.satellite_pos_z      = s.z;
        m.pseudorange_m        = rng + cb_truth;
        m.pseudorange_rate_m_s = 0.0; // static rx + static sats + zero clock drift
        m.elevation_rad        = look[i].el_deg * DEG;
        meas.push_back( m );
    }
    // One SV is a gross outlier (cross-correlation / false-preamble track): pseudorange off by 5000 km.
    meas[3].pseudorange_m += 5.0e6;

    Position_solver solver;
    const auto      out = solver.compute_solution( meas, 0.0 ); // seeds on the first call

    REQUIRE( out.has_value() ); // still seeds: RAIM drops the outlier, 6 clean SVs keep good geometry
    REQUIRE( out->valid );
    // Truth recovered to within a metre -> the 5000 km outlier was rejected (had it stayed in the LS the
    // seed would be hundreds of km off). The remaining 6 SVs are noise-free so the fit is essentially exact.
    REQUIRE( out->ecef_x_m == Catch::Approx( user.x ).margin( 1.0 ) );
    REQUIRE( out->ecef_y_m == Catch::Approx( user.y ).margin( 1.0 ) );
    REQUIRE( out->ecef_z_m == Catch::Approx( user.z ).margin( 1.0 ) );
    REQUIRE( out->clock_bias_m == Catch::Approx( cb_truth ).margin( 1.0 ) );
}
#endif
