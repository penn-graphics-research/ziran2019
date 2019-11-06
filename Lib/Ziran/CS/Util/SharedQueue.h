#ifndef SHARED_QUEUE_H
#define SHARED_QUEUE_H
#include <signal.h>
#include <condition_variable>
#include <exception>
#include <mutex>
#include <queue>
#include <iostream>

namespace ZIRAN {
/** Multiple producer, multiple consumer thread safe queue
 * Since 'return by reference' is used this queue won't throw */
template <typename T>
class SharedQueue {
    std::queue<T> queue;
    mutable std::mutex m;
    std::condition_variable data_cond;

    SharedQueue& operator=(const SharedQueue&) = delete;
    SharedQueue(const SharedQueue& other) = delete;

public:
    SharedQueue() {}

    void push(T&& item)
    {
        {
            std::lock_guard<std::mutex> lock(m);
            queue.emplace(std::move(item));
        }
        data_cond.notify_one();
    }

    /// \return immediately, with true if successful retrieval
    bool tryAndPop(T& popped_item)
    {
        std::lock_guard<std::mutex> lock(m);
        if (queue.empty()) {
            return false;
        }
        popped_item = std::move(queue.front());
        queue.pop();
        return true;
    }

    /// Try to retrieve, if no items, wait till an item is available and try again
    void waitAndPop(T& popped_item)
    {
        std::unique_lock<std::mutex> lock(m);
        while (queue.empty()) {
            data_cond.wait(lock);
            //  This 'while' loop is equal to
            // data_cond.wait(lock, [](bool result){return !queue.empty();});
        }
        popped_item = std::move(queue.front());
        queue.pop();
    }

    bool empty() const
    {
        std::lock_guard<std::mutex> lock(m);
        return queue.empty();
    }

    unsigned size() const
    {
        std::lock_guard<std::mutex> lock(m);
        return queue.size();
    }
};
} // namespace ZIRAN
#endif
