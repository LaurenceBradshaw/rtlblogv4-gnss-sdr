#include "logging.h"

namespace
{
struct Log_entry
{
    logging::Level                        level;
    std::string                           message;
    std::chrono::system_clock::time_point timestamp;
};

struct Logger_state
{
    std::queue<Log_entry>   log_queue_;
    std::mutex              queue_mutex_;
    std::condition_variable queue_cv_;
    std::thread             logger_thread_;
    bool                    running_ { false };
    bool                    stop_requested_ { false };

    ~Logger_state()
    {
        // Ensure logger thread is stopped and joined on process exit to avoid
        // terminating the application if tests or other code forgot to call
        // `logging::stop_log()`.
        {
            std::lock_guard<std::mutex> lock( queue_mutex_ );
            stop_requested_ = true;
        }
        queue_cv_.notify_all();
        if( logger_thread_.joinable() )
        {
            try
            {
                logger_thread_.join();
            }
            catch( ... )
            {
                // swallow exceptions during static teardown
            }
        }
    }
};

Logger_state& get_logger_state()
{
    static Logger_state state;
    return state;
}

struct FileBuffer
{
    std::string data;
    std::size_t max_size { 8192 };
};

static std::unordered_map<std::string, FileBuffer> g_file_buffers;
static std::mutex                                  g_file_buffers_mutex;

void buffered_dump( const std::string& file_name, const std::string& content )
{
    std::lock_guard<std::mutex> lock( g_file_buffers_mutex );

    auto& fb = g_file_buffers[file_name];
    fb.data += content + "\n";

    if( fb.data.size() >= fb.max_size )
    {
        const std::filesystem::path p = logging::get_log_dir() / file_name;
        std::ofstream               ofs( p, std::ios::app );
        if( ofs )
        {
            ofs << fb.data;
            ofs.close();
            fb.data.clear();
        }
    }
}

void flush_all_buffers()
{
    std::lock_guard<std::mutex> lock( g_file_buffers_mutex );

    for( auto& kv : g_file_buffers )
    {
        const std::string& file_name = kv.first;
        FileBuffer&        fb        = kv.second;

        if( fb.data.empty() )
            continue;

        const std::filesystem::path p = logging::get_log_dir() / file_name;
        std::ofstream               ofs( p, std::ios::app );
        if( ofs )
        {
            ofs << fb.data;
            ofs.close();
            fb.data.clear();
        }
    }
}

void print_log( const Log_entry& entry )
{
    static const std::string log_file = "console_log.txt";
    std::ostream&            stream   = ( entry.level == logging::Level::Error ) ? std::cerr : std::cout;

    const std::time_t time_t_val = std::chrono::system_clock::to_time_t( entry.timestamp );
    const std::tm     tm         = *std::localtime( &time_t_val );
    char              time_buf[9];
    std::strftime( time_buf, sizeof( time_buf ), "%H:%M:%S", &tm );

    const std::string content = fmt::format( "{} [{}]: {}", time_buf, logging::to_string( entry.level ), entry.message );

    stream << fmt::format(
        "{}{}{}\n",
        logging::to_string( logging::to_colour( entry.level ) ),
        content,
        logging::to_string( logging::Colour::Reset )
    );

    // buffered_dump( log_file, content );
}

void logger_loop()
{
    Logger_state& state = get_logger_state();

    for( ;; )
    {
        Log_entry entry;

        {
            std::unique_lock<std::mutex> lock( state.queue_mutex_ );
            state.queue_cv_.wait( lock, [&state] { return state.stop_requested_ || !state.log_queue_.empty(); } );

            if( state.stop_requested_ && state.log_queue_.empty() )
            {
                break;
            }

            entry = std::move( state.log_queue_.front() );
            state.log_queue_.pop();
        }

        print_log( entry );
    }
}

void ensure_logger_started()
{
    Logger_state&               state = get_logger_state();
    std::lock_guard<std::mutex> lock( state.queue_mutex_ );

    if( state.running_ )
    {
        return;
    }

    state.stop_requested_ = false;
    state.running_        = true;
    state.logger_thread_  = std::thread( logger_loop );
}

void enqueue_log( Log_entry entry )
{
    Logger_state& state = get_logger_state();
    {
        std::lock_guard<std::mutex> lock( state.queue_mutex_ );
        state.log_queue_.push( std::move( entry ) );
    }

    state.queue_cv_.notify_one();
}

} // namespace

namespace logging
{

std::filesystem::path get_log_dir()
{
    using Clock      = std::chrono::system_clock;
    using Time_point = Clock::time_point;

    static const Time_point      timestamp  = Clock::now();
    static const std::time_t     time_t_val = Clock::to_time_t( timestamp );
    static const std::tm         tm         = *std::localtime( &time_t_val );
    static char                  time_buf[21];
    static bool                  buf_set = false;
    static std::filesystem::path log_dir;
    static std::filesystem::path exe_dir = std::filesystem::canonical( "/proc/self/exe" ).parent_path();
    if( !buf_set )
    {
        std::strftime( time_buf, sizeof( time_buf ), "%Y_%m_%d_%H_%M_%S_", &tm );
        buf_set = true;
        static std::string pid( std::to_string( getpid() ) );
        static std::string timestamp_dir( time_buf + pid );
        log_dir = exe_dir / "logs" / timestamp_dir;
        std::filesystem::create_directories( log_dir );
    }

    return log_dir;
}

void log( logging::Level level, const std::string& str )
{
#ifdef NDEBUG
    // Discard debug level messages for release builds
    if( level == logging::Level::Debug )
    {
        return;
    }
#endif
    ensure_logger_started();

    Log_entry entry { level, std::move( str ), std::chrono::system_clock::now() };
    enqueue_log( std::move( entry ) );
}

void stop_log()
{
    Logger_state& state = get_logger_state();

    std::thread logger_thread;
    {
        std::lock_guard<std::mutex> lock( state.queue_mutex_ );

        if( !state.running_ )
        {
            return;
        }

        state.stop_requested_ = true;
        logger_thread         = std::move( state.logger_thread_ );
    }

    state.queue_cv_.notify_all();

    if( logger_thread.joinable() )
    {
        logger_thread.join();
    }

    // Ensure any buffered file writes are flushed to disk before shutdown.
    // flush_all_buffers();

    std::lock_guard<std::mutex> lock( state.queue_mutex_ );
    state.running_        = false;
    state.stop_requested_ = false;
}

} // namespace logging

#ifdef ENABLE_UNIT_TESTS
#include <catch2/catch_test_macros.hpp>

TEST_CASE( "logging_level_to_string_is_stable", "[logging]" )
{
    REQUIRE( logging::to_string( logging::Level::Info ) == "INFO" );
    REQUIRE( logging::to_string( logging::Level::Warning ) == "WARNING" );
    REQUIRE( logging::to_string( logging::Level::Error ) == "ERROR" );
    REQUIRE( logging::to_string( logging::Level::Debug ) == "DEBUG" );
}

// Dumping currently disabled.
// TEST_CASE( "test_dump_writes_file", "[logging]" )
// {
//     const std::string content = "hello_from_test";

//     // Use the public API; buffered writes are flushed on stop_log().
//     logging::log( logging::Level::Info, content );
//     logging::stop_log();

//     const std::filesystem::path p = logging::get_log_dir() / "console_log.txt";
//     REQUIRE( std::filesystem::exists( p ) );

//     std::ifstream ifs( p );
//     REQUIRE( ifs.is_open() );
//     bool        found = false;
//     std::string line;
//     while( std::getline( ifs, line ) )
//     {
//         if( line.find( content ) != std::string::npos )
//         {
//             found = true;
//             break;
//         }
//     }
//     REQUIRE( found );
// }

TEST_CASE( "test_stop_idempotent", "[logging]" )
{
    // start by logging a message (lazy-starts the logger)
    logging::log( logging::Level::Info, "test message for idempotent stop" );

    REQUIRE_NOTHROW( logging::stop_log() );
    // calling stop again must be safe
    REQUIRE_NOTHROW( logging::stop_log() );
}

TEST_CASE( "test_multiple_and_stop", "[logging]" )
{
    for( int i = 0; i < 10; ++i )
    {
        logging::log( logging::Level::Debug, "iter_" + std::to_string( i ) );
    }

    // stop_log should flush queued messages and return without throwing
    REQUIRE_NOTHROW( logging::stop_log() );
}
#endif