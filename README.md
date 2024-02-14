# Control Systems and Adaptive Strategies in Rust

Rust lang. to build adaptive control mechanisms, for embedded systems, to balance double, triple or more pendulum system.
Author: [Shashank S.K](https://github.com/CosmicBug), [Ting-Wei C.](https://github.com/tim840818)



Control systems are, as the name suggests, a system that generates sequences of input (i.e., controls) to stabilize or bring the input-fed system of interest into the desired configuration. The system of interest and expected noise are generally modeled systematically before the engineered control system is deployed. However, in more real-world scenarios, it would be rather impractical to model all the possible deviations and noise and incorporate those into the design of the control systems. Adaptive control systems are designed to locally model and respond to the unexpected external inputs and noise to reach or regain the desired configurations of the system of interest.

This project implements the task of balancing an inverted double pendulum as a toy model in Rust, with the purpose of contrasting control system methods and the corresponding implemented language performance. With a special emphasis on executing in an embedded systems environment, which generally demands a lean hardware utility for a fast system response.
