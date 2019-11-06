#ifndef ACTIVE_H
#define ACTIVE_H
#include <Ziran/CS/Util/SharedQueue.h>
#include <functional>
#include <memory>
#include <thread>

namespace ZIRAN {

/**
  Active object which enqueues callbacks to be run on another
  thread.  Safe for concurrent use from multiple threads.
  The callbacks will be run in the order recieved.
  */
class Active {
private:
    Active(std::function<bool()>&& shouldStop);
    // Construction ONLY through factory create();

    Active(const Active&) = delete;
    Active& operator=(const Active&) = delete;

    void run();

    SharedQueue<std::function<void()>> mq;
    std::function<bool()> shouldStop;
    std::thread thd;
    bool done;

public:
    virtual ~Active()
    {
        send([this] { done = true; });
        thd.join();
    }

    void send(std::function<void()>&& msg);

    unsigned size() const;

    /// Factory: safe construction of object before thread start
    // The callback stop is used to tell if the active thread should stop early without
    // executing all the callbacks in the queue
    static std::unique_ptr<Active> create(std::function<bool()>&& stop = []() { return false; });
};
} // namespace ZIRAN
#endif
