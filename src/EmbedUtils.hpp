//
// Created by Andrew Bailey on 2019-06-18.
//

#ifndef EMBED_FAST5_UTILS_H
#define EMBED_FAST5_UTILS_H

//Boost Libraries
#include <boost/filesystem.hpp>
#include <boost/coroutine2/all.hpp>
#include <boost/range/iterator_range.hpp>

// Standard Libraries
#include <exception>
#include <string>
#include <sstream>
#include <iostream>
#include <set>
#include <iomanip>
#include <mutex>

using namespace boost::filesystem;
using namespace boost::coroutines2;
using namespace std;

typedef coroutine<path> dir_coro;


namespace embed_utils{
  bool are_characters_in_string(string &characters, string& my_string);
  size_t get_file_size(const path &filename);
  int64_t string_to_int(std::string& str_int);
  path make_dir(path &output_path);
  bool compare_files(const path &p1, const path &p2);
  bool copyDir(const path& source, const path& destination);
  std::vector<std::string> split_string(string& in, char delimiter);
  float string_to_float(const string &str_int);
  string sort_string(string &str);
  vector<string> all_lexicographic_recur(string &characters, string &data, uint64_t last, uint64_t index);
  vector<string> all_string_permutations(const string &characters, uint64_t &length);
  string remove_duplicate_characters(const string &input_string);
  void dir_iterator_coroutine(dir_coro::push_type& yield, path& directory, string& ext);
  dir_coro::pull_type list_files_in_dir(path& directory, string& ext);
  std::map<string, string> create_ambig_bases();
  std::map<string, string> create_ambig_bases2(string config_file);
  tuple<uint64_t, uint64_t, uint64_t, uint64_t> get_time(std::function<void()> bound_function);
  string get_time_string(std::function<void()> bound_function);
  int64_t lines_in_file(path &file_path);
  uint64_t number_of_columns(const path &file_path, char sep='\t');
  std::set<char> add_string_to_set(const std::set<char>& a, const string& b);
  string char_set_to_string(std::set<char> a);
  std::set<char> string_to_char_set(const string& a);
  path make_dir(path &output_path);
  uint64_t compute_string_hash(string const& s);
  /**
  * Remove all empty file paths from vector
  *
  * @param file_paths: vector of paths to files
  * @param ext: extension to keep. ".tsv"
  */
  template<class T>
  vector<path> filter_emtpy_files(vector<T>& file_paths, string ext){
    vector<path> all_files;
    for (auto& a_string: file_paths) {
      path a_path(a_string);
  //        filter for files that are regular, end with tsv and are not empty
      if (is_regular_file(a_path) and a_path.extension().string() == ext and get_file_size(a_path.string()) > 0) {
        all_files.push_back(a_path);
      }
    }
    return all_files;
  }
  /**
  * Capture cout and cerr
  * source: https://stackoverflow.com/questions/5419356/redirect-stdout-stderr-to-a-string
  */
  class Redirect {

   public:
    bool out;
    bool err;

    Redirect(bool out=true, bool err=false): out(out), err(err) {
      if (out){
        coutold = std::cout.rdbuf( coutbuffer.rdbuf() ); // redirect cout to buffer stream
      }
      if (err){
        cerrold = std::cerr.rdbuf( cerrbuffer.rdbuf() ); // redirect cerr to buffer stream
      }

    }

    std::string get_cout() {
      return coutbuffer.str(); // get string
    }
    std::string get_cerr() {
      return cerrbuffer.str(); // get string
    }


    ~Redirect() {
      if (out){
        std::cout.rdbuf(coutold); // reverse redirect
      }
      if (err){
        std::cerr.rdbuf(cerrold); // reverse redirect
      }
    }

   private:
    std::stringstream coutbuffer;
    std::stringstream cerrbuffer;
    std::streambuf * coutold;
    std::streambuf * cerrold;
  };

//  https://codereview.stackexchange.com/questions/186535/progress-bar-in-c
  class ProgressBar {
    static const auto overhead = sizeof " [100%]";

    std::ostream &os;
    const std::size_t bar_width;
    std::string message;
    const std::string full_bar;
    std::mutex _mutex;

   public:
    ProgressBar(std::ostream &os, std::size_t line_width,
                std::string message_, const char symbol = '.')
        : os{os},
          bar_width{line_width - overhead},
          message{std::move(message_)},
          full_bar{std::string(bar_width, symbol) + std::string(bar_width, ' ')} {
      if (message.size() + 1 >= bar_width || message.find('\n') != message.npos) {
        os << message << '\n';
        message.clear();
      } else {
        message += ' ';
      }
      write(0.0);
    }

    // not copyable
    ProgressBar(const ProgressBar &) = delete;
    ProgressBar &operator=(const ProgressBar &) = delete;

    ~ProgressBar() {
      write(1.0);
      os << '\n';
    }

    void write(double fraction) {
      std::unique_lock<std::mutex> lock(_mutex);
      // clamp fraction to valid range [0,1]
      if (fraction < 0)
        fraction = 0;
      else if (fraction > 1)
        fraction = 1;

      auto width = bar_width - message.size();
      auto offset = bar_width - static_cast<unsigned>(width * fraction);

      os << '\r' << message;
      os.write(full_bar.data() + offset, width);
      os << " [" << std::setw(3) << static_cast<int>(100 * fraction) << "%] " << std::flush;
      _mutex.unlock();
    }
  };
}

/// Exception type for assertion failures
class AssertionFailureException : public std::exception
{
 private:
  const char* expression;
  const char* file;
  int line;
  std::string message;
  std::string report;

 public:

  /// Helper class for formatting assertion message
  class StreamFormatter
  {
   private:

    std::ostringstream stream;

   public:

    operator std::string() const
    {
      return stream.str();
    }

    template<typename T>
    StreamFormatter& operator << (const T& value)
    {
      stream << value;
      return *this;
    }
  };

  /// Log error before throwing
  void LogError()
  {
#ifdef THROWASSERT_LOGGER
    THROWASSERT_LOGGER(report);
#else
    std::cerr << report << std::endl;
#endif
  }

  /// Construct an assertion failure exception
  AssertionFailureException(const char* expression, const char* file, int line, const std::string& message)
      : expression(expression)
      , file(file)
      , line(line)
      , message(message)
  {
    std::ostringstream outputStream;

    if (!message.empty())
    {
      outputStream << message << ": ";
    }

    std::string expressionString = expression;

    if (expressionString == "false" || expressionString == "0" || expressionString == "FALSE")
    {
      outputStream << "Unreachable code assertion";
    }
    else
    {
      outputStream << "Assertion '" << expression << "'";
    }

    outputStream << " failed in file '" << file << "' line " << line;
    report = outputStream.str();

    LogError();
  }

  /// The assertion message
  const char* what() const noexcept override
  {
    return report.c_str();
  }

  /// The expression which was asserted to be true
  const char* Expression() const noexcept
  {
    return expression;
  }

  /// Source file
  const char* File() const noexcept
  {
    return file;
  }

  /// Source line
  int Line() const noexcept
  {
    return line;
  }

  /// Description of failure
  const char* Message() const noexcept
  {
    return message.c_str();
  }

  ~AssertionFailureException() noexcept override
  = default;
};

/// Assert that EXPRESSION evaluates to true, otherwise raise AssertionFailureException with associated MESSAGE (which may use C++ stream-style message formatting)
#define throw_assert(EXPRESSION, MESSAGE) if(!(EXPRESSION)) { throw AssertionFailureException(#EXPRESSION, __FILE__, __LINE__, (AssertionFailureException::StreamFormatter() << MESSAGE)); }

#endif //EMBED_FAST5_UTILS_H
