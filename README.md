# Multi-Commodity Flows over Time
## Requirements:
[Python](https://python.org/) >= 3.7

[numpy](https://numpy.org/) >= 1.17

[networkx](https://networkx.github.io/) >= 2.4

## Usage:
    usage: mainControl.py [-h] networkFile inflowFile
    Compute multi-commodity Flows over Time.
    
    positional arguments:
      networkFile  Spillback network in .cg-Format.
      inflowFile   Inflow .txt File containing the commodities.
    
    optional arguments:
      -h, --help   show this help message and exit

## Example:

    python mainControl.py example/simple_merge.cg example/inflow.txt 
    Starting computation.
    
    Iteration 1: Node 1 | Theta 1.00 | Alpha 1.00 | 4 nodes in queue
    Iteration 21: Node 2 | Theta 12.00 | Alpha 3.00 | 4 nodes in queue
    Iteration 41: Node t | Theta 37.00 | Alpha 2.00 | 2 nodes in queue
    Iteration 45: Node t | Theta 46.50 | Alpha inf | 0 nodes in queue
    Done after 45 iterations. Elapsed time: 0.02 seconds.
    Generating output.
    
    Finished generating output. Total elapsed: 0.02 seconds
