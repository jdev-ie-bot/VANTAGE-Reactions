#include "include/mock_particle_group.hpp"
#include "include/test_common.hpp"
#include <cmath>
#include <limits>

using namespace VANTAGE::Reactions;

template <typename T> class CommonStrategies : public ::testing::Test {
protected:
  NP::ParticleGroupSharedPtr group;
  NP::Sym<T> sym{"COMMON_VALUES"};

  void SetUp() override {
    group = create_test_particle_group(96);
    group->add_particle_dat(sym, 2);
    group->add_particle_dat(NP::Sym<INT>("COMMON_FLAG"), 1);
    NP::particle_loop(
        "initialise_common_values", group,
        [](auto index, auto values) {
          values[0] = 7;
          values[1] = index.get_local_linear_index() % 3 - 1;
        },
        NP::Access::read(NP::ParticleLoopIndex{}), NP::Access::write(sym))
        ->execute();
  }

  void TearDown() override {
    group->sycl_target->free();
    group->domain->mesh->free();
  }

  template <typename PREDICATE>
  void check_marker(std::shared_ptr<MarkingStrategy> marker, PREDICATE expected,
                    int component = 0) {
    const auto flag = NP::Sym<INT>("COMMON_FLAG");
    NP::particle_loop(
        "clear_flags", group, [](auto flags) { flags[0] = 0; },
        NP::Access::write(flag))
        ->execute();
    auto selected = marker->make_marker_subgroup(NP::particle_sub_group(group));
    NP::particle_loop(
        "mark_selected", selected, [](auto flags) { flags[0] = 1; },
        NP::Access::write(flag))
        ->execute();

    INT expected_count = 0;
    for (int cell = 0; cell < group->domain->mesh->get_cell_count(); cell++) {
      auto values = group->get_cell(sym, cell);
      auto flags = group->get_cell(flag, cell);
      for (int row = 0; row < values->nrow; row++) {
        const bool matches = expected(values->at(row, component));
        EXPECT_EQ(flags->at(row, 0), matches);
        expected_count += matches;
      }
    }
    EXPECT_EQ(selected->get_npart_local(), expected_count);
  }
};

using PropertyTypes = ::testing::Types<INT, REAL>;
TYPED_TEST_SUITE(CommonStrategies, PropertyTypes);

TYPED_TEST(CommonStrategies, NamedComparisons) {
  // Component zero is constant: exercise default arguments and empty/full sets.
  this->check_marker(make_equal_to_marking_strategy(this->sym, 7),
                     [](auto value) { return value == 7; });
  this->check_marker(make_less_than_marking_strategy(this->sym, 7),
                     [](auto value) { return value < 7; });
  this->check_marker(make_greater_than_marking_strategy(this->sym, 6),
                     [](auto value) { return value > 6; });

  // A zero threshold must not reset a nonzero component to the default.
  this->check_marker(
      make_equal_to_marking_strategy(this->sym, 0, 1),
      [](auto value) { return value == 0; }, 1);
  this->check_marker(
      make_less_than_marking_strategy(this->sym, 0, 1),
      [](auto value) { return value < 0; }, 1);
  this->check_marker(
      make_greater_than_marking_strategy(this->sym, 0, 1),
      [](auto value) { return value > 0; }, 1);
}

TYPED_TEST(CommonStrategies, OtherComparisons) {
  this->check_marker(
      make_comparison_marking_strategy(this->sym, 0,
                                       std::less_equal<TypeParam>{}, 1),
      [](auto value) { return value <= 0; }, 1);
  this->check_marker(
      make_comparison_marking_strategy(this->sym, 0,
                                       std::greater_equal<TypeParam>{}, 1),
      [](auto value) { return value >= 0; }, 1);
  this->check_marker(
      make_comparison_marking_strategy(this->sym, 0,
                                       std::not_equal_to<TypeParam>{}, 1),
      [](auto value) { return value != 0; }, 1);
}

TYPED_TEST(CommonStrategies, CapturesValuesAndSurvivesFactoryScope) {
  TypeParam threshold = 0;
  int component = 1;
  auto marker =
      make_less_than_marking_strategy(this->sym, threshold, component);
  auto assignment =
      make_set_value_transformation_strategy(this->sym, threshold, component);
  threshold = 100;
  component = 0;
  this->check_marker(marker, [](auto value) { return value < 0; }, 1);

  std::vector<std::vector<TypeParam>> original_values;
  for (int cell = 0; cell < this->group->domain->mesh->get_cell_count();
       cell++) {
    auto values = this->group->get_cell(this->sym, cell);
    original_values.emplace_back();
    for (int row = 0; row < values->nrow; row++) {
      original_values.back().push_back(values->at(row, 1));
    }
  }

  // Transform only a subgroup, keeping the other component and particles
  // intact.
  auto selected =
      marker->make_marker_subgroup(NP::particle_sub_group(this->group));
  assignment->transform(selected);
  for (int cell = 0; cell < this->group->domain->mesh->get_cell_count();
       cell++) {
    auto values = this->group->get_cell(this->sym, cell);
    auto flags = this->group->get_cell(NP::Sym<INT>("COMMON_FLAG"), cell);
    for (int row = 0; row < values->nrow; row++) {
      EXPECT_EQ(values->at(row, 0), 7);
      if (flags->at(row, 0)) {
        EXPECT_EQ(values->at(row, 1), 0);
      } else {
        EXPECT_EQ(values->at(row, 1), original_values[cell][row]);
      }
    }
  }
}

TYPED_TEST(CommonStrategies, ComposeMarkersAndRemove) {
  auto wrapper = TransformationWrapper(
      std::vector<std::shared_ptr<MarkingStrategy>>{
          make_greater_than_marking_strategy(this->sym, -1, 1),
          make_less_than_marking_strategy(this->sym, 1, 1)},
      make_transformation_strategy<SimpleRemovalTransformationStrategy>());
  const auto original_count = this->group->get_npart_local();
  wrapper.transform(this->group);
  EXPECT_EQ(this->group->get_npart_local(), 2 * original_count / 3);
  for (int cell = 0; cell < this->group->domain->mesh->get_cell_count();
       cell++) {
    auto values = this->group->get_cell(this->sym, cell);
    for (int row = 0; row < values->nrow; row++) {
      EXPECT_NE(values->at(row, 1), 0);
    }
  }
  // Reusing a wrapper with no matching particles must be harmless.
  wrapper.transform(this->group);
  EXPECT_EQ(this->group->get_npart_local(), 2 * original_count / 3);
}

TYPED_TEST(CommonStrategies, AssignmentDefaultAndNonzeroComponents) {
  auto wrapper = TransformationWrapper(
      std::vector<std::shared_ptr<MarkingStrategy>>{
          make_equal_to_marking_strategy(this->sym, 0, 1)},
      make_set_value_transformation_strategy(this->sym, 42));
  wrapper.transform(this->group);
  for (int cell = 0; cell < this->group->domain->mesh->get_cell_count();
       cell++) {
    auto values = this->group->get_cell(this->sym, cell);
    for (int row = 0; row < values->nrow; row++) {
      EXPECT_EQ(values->at(row, 0), values->at(row, 1) == 0 ? 42 : 7);
    }
  }

  auto assignment = make_set_value_transformation_strategy(this->sym, -3, 1);
  assignment->transform(NP::particle_sub_group(this->group));
  for (int cell = 0; cell < this->group->domain->mesh->get_cell_count();
       cell++) {
    auto values = this->group->get_cell(this->sym, cell);
    for (int row = 0; row < values->nrow; row++) {
      EXPECT_EQ(values->at(row, 1), -3);
      EXPECT_TRUE(values->at(row, 0) == 7 || values->at(row, 0) == 42);
    }
  }
}

TYPED_TEST(CommonStrategies, RejectNegativeComponents) {
  if (std::getenv("TEST_NESOASSERT") == nullptr) {
    GTEST_SKIP() << "Set TEST_NESOASSERT=ON to check assertion failures.";
  }
  EXPECT_THROW(make_equal_to_marking_strategy(this->sym, 0, -1),
               std::logic_error);
  EXPECT_THROW(make_set_value_transformation_strategy(this->sym, 0, -1),
               std::logic_error);
}

TEST(CommonStrategiesReal, ExactEqualityAndUnorderedComparisons) {
  const auto nan = std::numeric_limits<REAL>::quiet_NaN();
  // Exercise floating-point semantics on device, including NaN.
  auto group = create_test_particle_group(96);
  auto sym = NP::Sym<REAL>("WEIGHT");
  make_set_value_transformation_strategy(sym, nan)->transform(
      NP::particle_sub_group(group));
  for (auto marker : {make_equal_to_marking_strategy(sym, nan),
                      make_less_than_marking_strategy(sym, 0),
                      make_greater_than_marking_strategy(sym, 0)}) {
    EXPECT_EQ(marker->make_marker_subgroup(NP::particle_sub_group(group))
                  ->get_npart_local(),
              0);
  }
  make_set_value_transformation_strategy(sym, 1.0)->transform(
      NP::particle_sub_group(group));
  EXPECT_EQ(make_equal_to_marking_strategy(sym, std::nextafter(1.0, 2.0))
                ->make_marker_subgroup(NP::particle_sub_group(group))
                ->get_npart_local(),
            0);
  group->sycl_target->free();
  group->domain->mesh->free();
}
