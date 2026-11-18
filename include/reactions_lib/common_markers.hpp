#ifndef REACTIONS_COMMON_MARKERS_H
#define REACTIONS_COMMON_MARKERS_H
#include "particle_properties_map.hpp"
#include "reactions/neso_particles_namespace_alias.hpp"
#include "transformation_wrapper.hpp"
#include <functional>
#include <type_traits>

namespace VANTAGE::Reactions {

/**
 * @brief Compare one component of an INT or REAL particle property with a
 * constant using a device-copyable binary predicate.
 *
 * For example, pass std::less_equal<REAL>{} for an inclusive upper bound.
 * Values and the predicate are copied into the strategy. Floating-point
 * comparisons use ordinary C++ rules, with no implicit tolerance.
 *
 * @param sym Property to read.
 * @param value Constant to compare against (converted to the property's type).
 * @param comparison Predicate called as comparison(property[component], value).
 * @param component Component to read; must be a valid index in the property.
 * @return A marking strategy that can be composed in a TransformationWrapper.
 */
template <typename T, typename COMPARISON>
inline std::shared_ptr<MarkingStrategy>
make_comparison_marking_strategy(NP::Sym<T> sym, std::common_type_t<T> value,
                                 COMPARISON comparison, int component = 0) {
  static_assert(std::is_same_v<T, INT> || std::is_same_v<T, REAL>,
                "Comparison marking requires an INT or REAL property.");
  NESOASSERT(component >= 0, "Comparison component must be non-negative.");
  return make_direct_marking_strategy(
      "comparison_marker",
      [=](auto property) { return comparison(property[component], value); },
      NP::Access::read(sym));
}

/**
 * @brief Mark particles whose selected component equals value exactly.
 * @param sym Property to read.
 * @param value Constant to compare against.
 * @param component Valid component index, defaulting to zero.
 */
template <typename T>
inline std::shared_ptr<MarkingStrategy>
make_equal_to_marking_strategy(NP::Sym<T> sym, std::common_type_t<T> value,
                               int component = 0) {
  return make_comparison_marking_strategy(sym, value, std::equal_to<T>{},
                                          component);
}

/**
 * @brief Mark particles whose selected component is strictly less than value.
 * @param sym Property to read.
 * @param value Constant to compare against.
 * @param component Valid component index, defaulting to zero.
 */
template <typename T>
inline std::shared_ptr<MarkingStrategy>
make_less_than_marking_strategy(NP::Sym<T> sym, std::common_type_t<T> value,
                                int component = 0) {
  return make_comparison_marking_strategy(sym, value, std::less<T>{},
                                          component);
}

/**
 * @brief Mark particles whose selected component is strictly greater than
 * value.
 * @param sym Property to read.
 * @param value Constant to compare against.
 * @param component Valid component index, defaulting to zero.
 */
template <typename T>
inline std::shared_ptr<MarkingStrategy>
make_greater_than_marking_strategy(NP::Sym<T> sym, std::common_type_t<T> value,
                                   int component = 0) {
  return make_comparison_marking_strategy(sym, value, std::greater<T>{},
                                          component);
}

/**
 * @brief Marking strategy that selects only those particles in cells containing
 * some minimum number of particles
 *
 */
struct MinimumNPartInCellMarker : MarkingStrategy {

public:
  MinimumNPartInCellMarker() = delete;

  /**
   * @brief Constructor for MinimumNPartInCellMarker.
   *
   * @param min_npart Minimum number of particles in a cell.
   */
  MinimumNPartInCellMarker(INT min_npart);

  /**
   * @brief Helper function to construct a particle sub group that's generated
   * via application of a marking strategy that selects only those particles in
   * cells containing min_npart particles.
   *
   * @param particle_group The particle group to apply the marking strategy to.
   *
   * @return A shared pointer of a ParticleSubGroup that contains the selected
   * particles.
   */
  NP::ParticleSubGroupSharedPtr
  make_marker_subgroup_v(NP::ParticleSubGroupSharedPtr particle_group);

private:
  INT min_npart;
};

/**
 * @brief Marking strategy that selects only those particles with a panic flag >
 * 0
 *
 */
struct PanickedParticleMarker : MarkingStrategy {

public:
  PanickedParticleMarker() = delete;

  /**
   * @brief Constructor for PanickedParticleMarker.
   *
   * @param properties_map (Optional) A std::map<int, std::string> object to be
   * used to remap the NP::Sym for the Panic property.
   */
  PanickedParticleMarker(
      const std::map<int, std::string> &properties_map = get_default_map());

  /**
   * @brief Helper function to construct a particle sub group that's generated
   * via application of a marking strategy that selects only those particles
   * that have panicked.
   *
   * @param particle_group The particle group to apply the marking strategy to.
   *
   * @return A shared pointer of a ParticleSubGroup that contains the selected
   * particles.
   */
  NP::ParticleSubGroupSharedPtr
  make_marker_subgroup_v(NP::ParticleSubGroupSharedPtr particle_group);

private:
  NP::Sym<INT> panic_sym;
};

/**
 * @brief Helper function to check if there are any panicked particles in a
 * particle group.
 *
 * @param particle_group The particle group that contains the particles to be
 * counted.
 * @param properties_map (Optional) A std::map<int, std::string> object to be
 * used to remap the NP::Sym for the Panic property.
 *
 * @return Boolean to indicate if any panicked particles exist in the particle
 * group.
 */
bool panicked(
    NP::ParticleSubGroupSharedPtr particle_group,
    const std::map<int, std::string> &properties_map = get_default_map());

} // namespace VANTAGE::Reactions

#endif
