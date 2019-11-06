#include <Ziran/CS/Util/Active.h>

namespace ZIRAN {

Active::Active(std::function<bool()>&& shouldStop)
    : shouldStop(shouldStop)
    , done(false) // Construction ONLY through factory create();

{
}

void Active::run()
{
    while (!done) {
        if (shouldStop())
            break;
        std::function<void()> func;
        mq.waitAndPop(func);
        func();
    }
}

void Active::send(std::function<void()>&& msg)
{
    mq.push(std::move(msg));
}

unsigned Active::size() const
{
    return mq.size();
}

std::unique_ptr<Active>
Active::create(std::function<bool()>&& stop)
{
    std::unique_ptr<Active> aPtr(new Active(std::move(stop)));
    aPtr->thd = std::thread(&Active::run, aPtr.get());
    return aPtr;
}
} // namespace ZIRAN
