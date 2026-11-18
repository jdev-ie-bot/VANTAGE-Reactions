void common_strategies_example(NP::ParticleGroupSharedPtr particle_group) {
  const auto weight = NP::Sym<REAL>("WEIGHT");
  const auto state = NP::Sym<INT>("INTERNAL_STATE");

  // Select species 1 with a low weight, then remove those particles.
  auto remove_low_weight = TransformationWrapper(
      std::vector<std::shared_ptr<MarkingStrategy>>{
          make_equal_to_marking_strategy(state, 1),
          make_less_than_marking_strategy(weight, 1e-12)},
      make_transformation_strategy<SimpleRemovalTransformationStrategy>());
  remove_low_weight.transform(particle_group);

  // Other comparisons can use standard function objects without a custom
  // lambda.
  auto non_positive =
      make_comparison_marking_strategy(weight, 0, std::less_equal<REAL>{});

  // Set a value only on the selected subgroup.
  auto reset_weight = TransformationWrapper(
      std::vector<std::shared_ptr<MarkingStrategy>>{non_positive},
      make_set_value_transformation_strategy(weight, 0));
  reset_weight.transform(particle_group);

  // An optional final component index selects vector-valued property entries.
  auto reset_velocity_y =
      make_set_value_transformation_strategy(NP::Sym<REAL>("VELOCITY"), 0, 1);
  reset_velocity_y->transform(NP::particle_sub_group(particle_group));
}
