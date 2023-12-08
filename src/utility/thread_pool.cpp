#include "thread_pool.h"

// initialize static variables
size_t ThreadPool::_numThreads = 0;
std::vector<ThreadPool::WorkerThread> ThreadPool::_workers;
ThreadPool::WorkerQueue ThreadPool::_tasks;
std::mutex ThreadPool::_terminateMutex;
bool ThreadPool::_terminate;


// function definitions
void ThreadPool::Go(size_t thread_id) {
    while (true) {
        std::function<void(int)> task;
        std::promise<void> promise;
        {
            std::unique_lock<std::mutex> lock(_tasks._mutex);

            _tasks._condition.wait(lock, []() {
                return !_tasks.empty() || _terminate;				// if we woke-up by accident check: are there tasks or should we quit
                });

            if (_tasks.empty()) return;								// dead

            task = _tasks.front().first;
            promise = std::move(_tasks.front().second);
            _tasks.pop();
        }

        task(thread_id);											// so do it
        promise.set_value();										// done!
    }
}
bool ThreadPool::Startup(size_t numThreads) {
    _terminate = false;
    _numThreads = numThreads;
    _workers.resize(_numThreads);
    for (unsigned i = 0; i < _numThreads; i++) {
        _workers[i]._id = i;								// set id
        _workers[i]._thread = std::thread(ThreadPool::Go, i);		// main loop
    }

    return true;
}
void ThreadPool::Shutdown() {
    {
        std::unique_lock<std::mutex> lock(_terminateMutex);
        _terminate = true; // use this flag in condition.wait
    }

    _tasks._condition.notify_all(); // wake up all threads.

    // Join all threads.
    for (auto& w : _workers)
        w._thread.join();

    _workers.clear();
}
std::future<void> ThreadPool::PushTask(std::function<void(int)> task)
{
    std::future<void> f;						// return value

    {
        std::unique_lock<std::mutex> lock(_tasks._mutex);
        _tasks.push(std::make_pair(task, std::move(std::promise<void>())));
        f = _tasks.back().second.get_future();
    }

    _tasks._condition.notify_one();

    return f;
}
/*
void ThreadPool::WaitTillFinished() {
{
}
*/