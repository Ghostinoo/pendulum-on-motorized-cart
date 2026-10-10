# Pendulum on a Motorized Cart

A small Pygame + SciPy simulator for a pendulum attached to a motorized cart. It renders the cart, wheels, and rod, while integrating the nonlinear equations of motion with `scipy.integrate.ode` using the `dopri5` solver. It lets the user control the disturb and easily write a controller to test it out.

## Requirements

- Python **3.12+**
- [`pygame`](https://www.pygame.org/)
- [`scipy`](https://scipy.org/)

## Installation
Create a folder for this project and move the files `requirements.txt` and `pendulum.py` into it. Then run the following commands:
### Windows
```cmd
python -m venv venv
.\venv\Scripts\activate
pip install -r requirements.txt
```
### Linux
```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Quick Start
Create a python file, for example main.py, in the same folder.
```py
# Import the class Pendulum
from pendulum import Pendulum

# Define the parameters needed to run the simulation.
PARAMS = {
    'R': 0.05, # radius of the wheels (m)
    'L': 0.5, # length of the rod (m)
    'g': 9.81, # gravitational acceleration (m/s^2)
    'mC': 5, # the mass of the cart (Kg)
    'mR': 2, # the mass of the rod (Kg)
    'mMW': 0.2, # the mass of the motorized wheel rotating assembly (Kg)
    'mW': 0.2, # the mass of the non-motorized wheel rotating assembly (Kg)
    'JR': 2 * ((2*0.5)**2) / 12, # the moment of inertia of the rod (taken from the COM) [Kg*m^2]
    'JW': 0.2 * (0.05**2) / 2, # the moment of inertia of the non-motorized wheel assembly (Kg*m^2)
    'JMW': 0.2 * (0.05**2) / 2
}

# Define the controller, if you wish.
def controller(state):
    x, dx, p, dp = state
    # Simple PD controller on pendulum angle and cart position.
    # Returns motor torque. NB: torque is positive anti-clockwise,
    # so negative feedback is achieved with positive coefficients.
    return 80*p + 20*dp + 1*x + 5*dx


sim = Pendulum(
    parameters=PARAMS,
    initial_state=(0.0, 0.0, 0.2, 0.0),  # x, dx, p, dp
)

sim.set_controller(controller)
sim.run()
```

## Controls

Key | Action
- | -
`ESC` | Close the simulation
`A` | Push the rod to the **left** from its end
`D` | Push the rod to the **right** from its end

## API
In order to understand the types you see in the function arguments, please go to `Types`.
### ```Pendulum(parameters, initial_state=None, initial_time=0, env_configuration=None)```
This constructor initializes the pendulum simulator.

- `parameters: TParameters` A python dictionary containing the system parameters. 

- `initial_state: TState` List or tuple containing the initial state of the system.

- `initial_time: float` Time at which you want the simulation to start. This must be **non-negative**.

- `env_configuration: TEnvironmentConfigurationPartial` A python dictionary containing the environment configuration.

---
### Properties
You can access these variables by dot notation on your class object. You are **not** supposed to write these. In order to change one of these, please re-call the constructor.
- `state: TState` The current state of the system.
- `parameters: TParameters` The current parameters of the system.
- `env_configuration: TEnvironmentConfiguration` Current environment configuration.

---
### Methods
These are the public methods that the class `Pendulum` supports:

- `set_controller(controller: TController) -> None` Sets the controller for the system. The `controller` function is called at every simulation step and it gets passed the state of the system. It is supposed to compute and **return as a `float`** the torque applied to one of the wheels.

- `run() -> None` Runs the real-time interactive simulation.
---
### Types

- **`TParameters`**

```py
type TParameters = {
    'R': float, # radius of the wheels [m]
    'L': float, # length of the rod [m]
    'g': float, # gravitational acceleration [m/s^2]
    'mC': float, # the mass of the cart [Kg]
    'mR': float, # the mass of the rod [Kg]
    'mMW': float, # the mass of the motorized wheel rotating assembly [Kg]
    'mW': float, # the mass of the non-motorized wheel rotating assembly [Kg]
    'JR': float, # the moment of inertia of the rod (taken from the COM) [Kg*m^2]
    'JW': float, # the moment of inertia of the non-motorized wheel assembly [Kg*m^2]
    'JMW': float # the moment of inertia of the motorized wheel assembly [Kg*m^2]
}
```

- **`TEnvironmentConfiguration`**

Note: if you see `TEnvironmentConfigurationPartial` instead, that just means you can write just a few attributes of this dict: there's no need to define every one of it.
```py
type TEnvironmentConfiguration = {
    'SCREEN_WIDTH': int # the width of the screen in pixels
    'SCREEN_HEIGHT': int # the height of the screen in pixels
    'FPS': int # the frames per second. Simulation interval = 1/FPS
    'BG_COLOR': tuple[int, int, int] # background color in HEX
    'FLOOR_COLOR': tuple[int, int, int] # floor color in HEX
    'CART_COLOR': tuple[int, int, int] # cart color in HEX
    'FLOOR_Y': int # floor height in pixels
    'CART_WIDTH_M': float # cart width (in meters)
    'CART_HEIGHT_M': float # cart height (in meters)
    'PX_TO_M': float # how many pixels correspond to one meter
    'DISTURBANCE_FORCE': float # the force applied to the end of the rod in Newtons when A or D are pressed.
}
```

- **`TState`**
```py
type TState = tuple(float, float, float, float)
# TState[0] = x  ,  TState[1] = dx
# TState[2] = p  ,  TState[3] = dp
```

- **`TController`**
```py
type TState = Callable[[TState], float]
# It is a function that accepts 1 argument of type TState and returns a float
```

## Final Notes
- The controller is not necessary, if you want to just simulate the system, just don't call `Pendulum.set_controller()`.
- `x` is positive when going right-to-left, `p` is positive when going anti-clockwise.
- The torque is applied to one of the two wheels, not to both (in order to simulate a real cart with one motor). If the torque is positive the wheel gets turned anti-clockwise.
- Control inputs and disturbance forces are assumed constant during each integration step. This means that if the simulation runs at `60 FPS`, then the period is `1/60`, thus the inputs are updated once every `1/60` of second. Change the `FPS` config option to reduce or enlarge this period. Higher FPS results in lower difference between the real-time controller implemented and the discrete one actually running.
- If integration fails, there's probably something wrong with your controller.
- The solver, in order to compute the state of the system, integrates numerically the differential equations that describe the system. It uses the `dopri5` algorithm with `rtol=1e-8` and `atol=1e-10`.