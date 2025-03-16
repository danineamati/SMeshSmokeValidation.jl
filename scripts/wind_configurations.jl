# Include all the disturbance settings in one place to avoid duplication

using SMeshSmokeValidation


DEFAULT_DISTURBANCE(Q) = FullDisturbances(
    HeightDistribution(1000.0 + 0.0 * Q, 80.0),
    WindDistribution(
        [0.0, -120, -180],
        [9.0, 5.0, 5.0],
        [0.7, 0.15, 0.15],
        5.0, 1.0;
        in_deg=true
    )
)


# https://mesonet.agron.iastate.edu/onsite/windrose/CA_ASOS/E16/E16_sep.png
HENRY_COE_NOMINAL(Q) = FullDisturbances(
    HeightDistribution(Q, 30.0),
    WindDistribution(
        [120.0, -50.0, -110.0],
        [12.0, 12.0, 2.0],
        [0.6, 0.3, 0.1],
        4.5, 3.0;
        in_deg=true
    )
)

HENRY_COE_FUZZED(Q) = FullDisturbances(
    HeightDistribution(Q*2, 5.0),
    WindDistribution(
        [120.0, -50.0, -110.0],
        [15.0, 9.0, 9.0],
        [0.15, 0.7, 0.15],
        5.5, 2.0;
        in_deg=true
    )
)

# This is the original configuration
# MALIBU_NOMINAL(Q) = FullDisturbances(
#     HeightDistribution(Q, 30.0),
#     WindDistribution(
#         [-120.0, 90.0, 150.0],
#         [9.0, 5.0, 5.0],
#         [0.1, 0.7, 0.2],
#         10.0, 5.0;
#         in_deg=true
#     )
# )
# https://mesonet.agron.iastate.edu/onsite/windrose/CA_ASOS/NTD/NTD_nov.png
MALIBU_NOMINAL(Q) = FullDisturbances(
    HeightDistribution(Q, 30.0),
    WindDistribution(
        [170.0, 40.0, 70.0],
        [7.0, 3.0, 5.0],
        [0.2, 0.6, 0.2],
        4.5, 3.0;
        in_deg=true
    )
)
# MALIBU_FUZZED(Q) = FullDisturbances(
#     HeightDistribution(Q*2, 5.0),
#     WindDistribution(
#         [-120.0, 90.0, 150.0],
#         [9.0, 5.0, 5.0],
#         [0.7, 0.15, 0.15],
#         5.0, 1.0;
#         in_deg=true
#     )
# )
MALIBU_FUZZED(Q) = FullDisturbances(
    HeightDistribution(Q*2, 5.0),
    WindDistribution(
        [170.0, 40.0, 70.0],
        [15.0, 9.0, 9.0],
        [0.7, 0.15, 0.15],
        5.5, 2.0;
        in_deg=true
    )
)

# Based on data from:
# https://mesonet.agron.iastate.edu/onsite/windrose/CA_ASOS/MHS/MHS_oct.png
SHASTA_NOMINAL(Q) = FullDisturbances(
    HeightDistribution(Q, 30.0),
    WindDistribution(
        [40.0, 120.0, -50.0],
        [9.0, 5.0, 2.0],
        [0.35, 0.4, 0.25],
        4.5, 3.0;
        in_deg=true
    )
)

SHASTA_FUZZED(Q) = FullDisturbances(
    HeightDistribution(Q*2, 5.0),
    WindDistribution(
        [40.0, 120.0, -50.0],
        [15.0, 9.0, 9.0],
        [0.15, 0.7, 0.15],
        5.5, 2.0;
        in_deg=true
    )
)



