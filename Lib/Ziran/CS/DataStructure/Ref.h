#ifndef REF_H
#define REF_H
namespace ZIRAN {
template <class T>
class Ref {
public:
    using Type = T;
    Ref(T& ref) noexcept
        : ptr(&ref) {}
    Ref(T&& ref) = delete;
    Ref(const Ref& ref) noexcept = default;

    Ref& operator=(const Ref& x) noexcept = default;

    operator T&() const noexcept { return *ptr; }
    T& get() const noexcept { return *ptr; }

    void set(const T& v) noexcept { *ptr = v; }

private:
    T* ptr;
};
template <class T>
Ref<T> makeRef(T& ref) noexcept
{
    return Ref<T>(ref);
}
template <class T>
Ref<T> makeRef(Ref<T> ref) noexcept
{
    return Ref<T>(ref);
}
template <class T>
void makeRef(const T&& ref) = delete;
template <class T>
Ref<const T> makeCRef(const T& ref) noexcept
{
    return Ref<const T>(ref);
}
template <class T>
Ref<const T> makeCRef(Ref<T> ref) noexcept
{
    return Ref<const T>(ref);
}
template <class T>
void makeCRef(const T&& ref) = delete;
} // namespace ZIRAN
#endif
