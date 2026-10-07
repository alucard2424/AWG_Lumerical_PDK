# AWG

Python tools for designing and routing **Arrayed Waveguide Gratings (AWGs)** in Lumerical.

## Idea


The main goal is to generate the AWG array while maintaining consistent geometry and waveguide spacing.

## Fixed-Bend Routing

The routing uses a **fixed bend radius** for all AWG channels.

Each route consists of:

* Input straight section
* Fixed-radius input bend
* Middle straight section
* Fixed-radius output bend
* Output straight section

The straight-section lengths are solved to achieve the required total optical path length for each channel.

For an AWG with path-length increment `delta_L`:

```text
L_i = L_0 + i × delta_L
```

where `L_0` is the reference channel length.

The resulting routes can then be drawn directly in Lumerical.
