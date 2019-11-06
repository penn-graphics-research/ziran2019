#ifndef PARTICLE_H
#define PARTICLE_H
#include <Ziran/CS/DataStructure/DataManager.h>

namespace ZIRAN {

template <class T, int dim>
class Particles : public DataManager {
public:
    using TV = Vector<T, dim>;

    // mass, X, V are assumed to have values for all entry ids
    DataArray<T>& mass;
    DataArray<TV>& X;
    DataArray<TV>& V;

    template <typename... Types>
    Particles(const AttributeName<Types>&... a)
        : DataManager(mass_name(), X_name(), V_name(), a...)
        , mass(get(mass_name()))
        , X(get(X_name()))
        , V(get(V_name()))
    {
    }

    virtual ~Particles() {}

    template <typename... DataTypes>
    DATA::DataAppender<T, TV, TV, DataTypes...> appender(const AttributeName<DataTypes>&... names)
    {
        return DataManager::appender(mass_name(), X_name(), V_name(), names...);
    }
    inline static AttributeName<T> mass_name()
    {
        return AttributeName<T>("m");
    }
    inline static AttributeName<TV> X_name()
    {
        return AttributeName<TV>("P");
    }
    inline static AttributeName<TV> V_name()
    {
        return AttributeName<TV>("V");
    }
};
} // namespace ZIRAN
#endif
