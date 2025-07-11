add_test([=[Autocalibration.Simple]=]  /Data/gautt/bundle_adjustment/autocalibration_project/build/tests/test_autocalibration [==[--gtest_filter=Autocalibration.Simple]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[Autocalibration.Simple]=]  PROPERTIES WORKING_DIRECTORY /Data/gautt/bundle_adjustment/autocalibration_project/build/tests SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
set(  test_autocalibration_TESTS Autocalibration.Simple)
