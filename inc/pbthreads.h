

#ifndef __PB_THREADS_H__
#define __PB_THREADS_H__
#include <thread>
#include <pthread.h>

class PBThread {
private:
  std::thread mCppThread;
  std::string mName;

public:
  PBThread() {mName = ""; }
  template<class Fn, class... Args>
  PBThread(std::string name, Fn&& fn, Args&&... args) : 
    mCppThread (std::thread(fn, args... )),
    mName(name) {
      pthread_setname_np(mCppThread.native_handle(), mName.c_str());
  }
  PBThread& operator=(PBThread&& newPBT){ 
    mName = newPBT.getThreadName(); 
    mCppThread = std::move(newPBT.getCppThread());
    return *this;
  }
  std::thread& getCppThread() { return mCppThread; }
  void join() { mCppThread.join();}
  bool joinable() { return mCppThread.joinable(); }
  std::string getThreadName() { return mName; }
};
#endif
