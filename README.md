# Multi-Commodity Flows Over Time
Given network and commodity information the Multi-Commodity Flow Over Time Tool computes a feasible flow over time as described in [[1]](#references).
This tool does not simulate any route-choice, as each commodity sends all its flow along a fixed path.
The output flow over time is described by the cumulative inflow functions for each commodity. Additionally, the travel time functions are provided for each commodity.

For creating the network the [NashFlowComputation tool](https://github.com/zimmer-m/NashFlowComputation) is needed.

## Requirements
[Python](https://python.org/) >= 3.7
## Usage
    usage: mainControl.py [-h] networkFile inflowFile
    Compute multi-commodity Flows over Time.
    
    positional arguments:
      networkFile  Network in .cg-Format.
      inflowFile   Inflow .txt File containing the commodities.
    
    optional arguments:
      -h, --help   show this help message and exit

## Example
An example is provided in [/example](/example). The files simple_merge.cg and inflow.txt are input files,
other files are the corresponding output files.

    python mainControl.py example/simple_merge.cg example/inflow.txt 
    Starting computation.
    
    Iteration 1: Node 1 | Theta 1.00 | Alpha 1.00 | 4 nodes in queue
    Iteration 21: Node 2 | Theta 12.00 | Alpha 3.00 | 4 nodes in queue
    Iteration 41: Node t | Theta 37.00 | Alpha 2.00 | 2 nodes in queue
    Iteration 45: Node t | Theta 46.50 | Alpha inf | 0 nodes in queue
    Done after 45 iterations. Elapsed time: 0.02 seconds.
    Generating output.
    
    Finished generating output. Total elapsed: 0.02 seconds

## References
[1] Theresa Ziemke, Leon Sering, Laura Vargas Koch, Max Zimmer, Kai Nagel and Martin Skutella. "Flows Over Time as Continuous Limits of Packet-Based Network Simulations." In: 23rd EURO Working Group on Transportation Meeting (EWGT), 2020. DOI: [TBD](https://svn.vsp.tu-berlin.de/repos/public-svn/publications/vspwp/2020/20-10/ZiemkeEtAl2020FlowsOverTimeAsLimitOfMATSim.pdf) 