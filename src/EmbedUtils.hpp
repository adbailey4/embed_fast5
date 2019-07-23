//
// Created by Andrew Bailey on 2019-06-18.
//

#ifndef EMBED_FAST5_UTILS_H
#define EMBED_FAST5_UTILS_H

#include <exception>
#include <string>
#include <sstream>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/coroutine2/all.hpp>

#include <boost/range/iterator_range.hpp>
#include <string>

using namespace boost::filesystem;
using namespace boost::coroutines2;
using namespace std;

typedef coroutine<path> dir_coro;


namespace embed_utils{
  bool are_characters_in_string(string &characters, string &my_string);
  size_t getFilesize(const std::string& filename);
  int64_t convert_to_int(std::string& str_int);
  path make_dir(path output_path);
  bool compareFiles(const std::string& p1, const std::string& p2);
  bool copyDir(path const & source, path const & destination);
  std::vector<std::string> split_string(string& in, char delimiter);
  std::vector<std::string> split_string2(string s, string delimiter, uint64_t size=16);
  float convert_to_float(std::string& str_int);
  string sort_string(string &str);
  vector<string> all_lexicographic_recur(string characters, string data, int last, int index);
  vector<string> all_string_permutations(string characters, int length);
  string remove_duplicate_characters(string input_string);
  void dir_iterator_coroutine(dir_coro::push_type& yield, path& directory, string& ext);
  dir_coro::pull_type list_files_in_dir(path& directory, string& ext);
  std::map<string, string> create_ambig_bases();
  tuple<uint64_t, uint64_t, uint64_t, uint64_t> get_time(std::function<void()> bound_function);
  string get_time_string(std::function<void()> bound_function);

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
  virtual const char* what() const throw()
  {
    return report.c_str();
  }

  /// The expression which was asserted to be true
  const char* Expression() const throw()
  {
    return expression;
  }

  /// Source file
  const char* File() const throw()
  {
    return file;
  }

  /// Source line
  int Line() const throw()
  {
    return line;
  }

  /// Description of failure
  const char* Message() const throw()
  {
    return message.c_str();
  }

  ~AssertionFailureException() throw()
  {
  }
};


/// Assert that EXPRESSION evaluates to true, otherwise raise AssertionFailureException with associated MESSAGE (which may use C++ stream-style message formatting)
#define throw_assert(EXPRESSION, MESSAGE) if(!(EXPRESSION)) { throw AssertionFailureException(#EXPRESSION, __FILE__, __LINE__, (AssertionFailureException::StreamFormatter() << MESSAGE)); }



#endif //EMBED_FAST5_UTILS_H
