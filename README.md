# Primary Vertex Reconstruction
## Analysis Tools

We will need different measures to verify how well the reconstruction is working.

- **Find the position of the “interesting” vertex (Signal event vertex):**
  - **Efficiency:** How many of the simulated signal vertices are reconstructed?
  - **Precision:** How close are the positions of the simulated and reconstructed vertices for signal events?
  - How often is the signal event also reconstructed as the signal event?
  
- **Identify other pile-up vertices:**
  - **Efficiency:** How many of the simulated vertices are reconstructed?
  - **Purity:** How many of the reconstructed vertices are incorrect? This controls for false positives and keeps the error rate small.
  - **Precision:** How close are the positions of the simulated and reconstructed vertices?
  - Where on the z-axis are the reconstructed and simulated vertices? What is the **distribution**?
  
- **Determine the number of pile-up vertices event-by-event:** How many events are there for each simulation?

- **Track--Vertex assignment:** For each simulated and reconstructed pair, how many of their associated tracks match (in percentage)?

## More Analysis Stuff, but for future use

Due to us wanting to shirten time it takes to reconstruct verteces, we will need to add some tools that help us verify how much time it takes to run the DA.
By having some graphs which show time versus number uf verteces or time versus number of steps in DA.

Additionally, we should also focus on comparing the "segmented" DA to "original" DA for vertex reconstruction. Here it could be useful to run the same analytics for both and compare them in one plot.


