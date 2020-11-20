# Tests

The following are integration tests run on CI at every commit but can also be run locally.

## How to run tests

Make sure you have the required Python libraries installed on your system (see [DEVELOP.md](../DEVELOP.md)), then from you command line, run the following command:

```
python3 run_tests.py <branch_a> <branch_b> <build_type>
```
Where `<branch_a>` and `<branch_b>` are the two branches you want to compare and `<build_type>` is either `Debug` or `Release`.

E.g.
```
pip3 install -r requirements.txt
python3 run_tests.py dmey/external-libs dmey/comp-tests Release
```

Tests outputs are saved under `tests/outputs` and currently include a boxplot of approximate errors for `u`, `v`, and `w`.

These tests are already automatically executed at each commit through CI.
