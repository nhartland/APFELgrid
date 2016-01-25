// The MIT License (MIT)

// Copyright (c) 2016 Stefano Carrazza, Luigi Del Debbio, Nathan Hartland

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <exception>
#include <stdexcept>
#include <string>

namespace NNPDF
{
  /**
   *
   * Simple and generic exception handlers
   *  - implements runtime and logic errors
   *
   */

  /// Error to be thrown when runtime error is detected
  class RuntimeException: public std::runtime_error
  {
  public:
    RuntimeException(const std::string& tag, const std::string& what) : std::runtime_error(std::string("[") + std::string(tag) + std::string("] error: ") + std::string(what)) {}
  };

  /// Error to be thrown when logic error is detected
  class LogicException: public std::logic_error
  {
  public:
    LogicException(const std::string& tag, const std::string& what) : std::logic_error(std::string("[") + std::string(tag) + std::string("] error: ") + std::string(what)) {}
  };

  //_______________________________________________________________________
  class FileError: public RuntimeException
  {
  public:
    FileError(const std::string& tag, const std::string& what) : RuntimeException(tag,what) {}
  };

  class EvaluationError: public RuntimeException
  {
  public:
    EvaluationError(const std::string& tag, const std::string& what) : RuntimeException(tag,what) {}
  };

  class InitError: public RuntimeException
  {
  public:
    InitError(const std::string& tag, const std::string& what) : RuntimeException(tag,what) {}
  };

  class RangeError: public LogicException
  {
  public:
    RangeError(const std::string& tag, const std::string& what) : LogicException(tag,what) {}
  };

  class LengthError: public LogicException
  {
  public:
    LengthError(const std::string& tag, const std::string& what) : LogicException(tag,what) {}
  };

  class LogError: public LogicException
  {
  public:
    LogError(const std::string& tag, const std::string& what) : LogicException(tag,what) {}
  };

  class UserError: public LogicException
  {
  public:
    UserError(const std::string& tag, const std::string& what) : LogicException(tag,what) {}
  };

}
