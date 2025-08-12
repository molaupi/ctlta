#pragma once

#include <chrono>
#include <cstdint>

// A timer to measure how long some code takes to execute.
class Timer {
 public:
  // Constructs a timer and starts it.
  Timer() : startTime(std::chrono::steady_clock::now()) {}

  // Returns the time elapsed since the timer was started.
  template <typename UnitT = std::chrono::milliseconds>
  int64_t elapsed() const {
    const auto now = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<UnitT>(now - startTime).count();
  }

  // Restarts the timer.
  void restart() {
    startTime = std::chrono::steady_clock::now();
  }

 private:
  std::chrono::steady_clock::time_point startTime; // Time point when the timer was started.
};


class NoOpTimer {
    public:
    // Constructs a no-op timer.
    NoOpTimer() {}

    // Returns the time elapsed since the timer was started, which is always zero.
    template <typename = std::chrono::milliseconds>
    int64_t elapsed() const {
        return 0;
    }

    // Restarts the timer, which does nothing in this case.
    void restart() {}
};