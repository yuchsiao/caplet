To install with default settings, 

1. Execute `install.sh`. Make sure `qmake` is in your system path.
2. Put a symbolic link or binary of FASTCAP in the same folder of `caplet_geo` for computing reference capacitance matrices.

To modify settings, see the file `caplet_solver/include/caplet_parameter.h` for `caplet_solver`.

To uninstall, execute `uninstall.sh`.

If you encounter any problem, go to the individual folder and then `make`.
If problems happen at `caplet_geo`, run `qmake` first and then `make`.
