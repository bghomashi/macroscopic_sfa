#ifndef __THREADPOOL_H__
#define __THREADPOOL_H__

#include <thread>
#include <vector>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <functional>
#include <future>
#include <complex>


class ThreadPool {
    typedef std::pair<std::function<void(size_t)>, std::promise<void>> WorkerTask;

    struct WorkerQueue : public std::queue<WorkerTask> {
        std::mutex _mutex;
        std::condition_variable _condition;
    };
    struct WorkerThread {
        std::thread _thread;
        size_t _id;
    };

    static size_t _numThreads;
    static std::vector<WorkerThread> _workers;
    static WorkerQueue _tasks;
    static std::mutex _terminateMutex;
    static bool _terminate;

    static void Go(size_t thread_id);
public:
    inline static size_t WorkerCount() {
        return _numThreads;
    }

    static bool Startup(size_t numThreads = std::thread::hardware_concurrency());
    static void Shutdown();

    static std::future<void> PushTask(std::function<void(int)> task);
    static void WaitTillFinished();
};


#endif