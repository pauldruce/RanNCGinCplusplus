# TODO for RanNCG Library

- [x] Make the `Action()` method in Simulation.h/c to be virtual so users can override this function with a more optimised calculation of the action.
- [x] Add some default override methods for `Action()` for type (1,1) geometries etc.
- [ ] Add some tests to check the optimised action calculation matches the default behaviour.
- [ ] Need to create a logging mechanism for errors. Currently some values of g2 seem to cause the simulation to make the step size very small and the acceptance_rate also becomes 0. Need to isolate these issues and log them for further investigation. So we don't waste compute time.
