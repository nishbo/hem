Neuron types:
0: Leaky i-a-f
1: Hodgkin-Huxley

Synapse types:
0: Static
1: Tsodyks-Markram
2: STDP with g-transmission function
3: STDP with Tsodyks-Markram
4: Tsodyks-Markram with g-transmission function as seen in NEST
5: Canonical TM with exc synapse STDP and no 2 multiplier
51: G First type
52: G Second type
53: G First type with cut
54: G Second type with cut

Topology types:
0: random network
1: sm-w network
2: m synapses go from neuron number 0

Delay types:
0: random (0, delay_max] msec
1: distance-dependent
2: minimal

Stimulation types:
0: Uniform initial current distribution (Imin; Imax)
1: Normal initial current distribution N(Imean, Isd) from Imin to Imax
2: As in 1, but if I<15pA, I=0
