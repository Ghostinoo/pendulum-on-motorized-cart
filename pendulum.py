from typing import Final, TypedDict, Callable
from scipy.integrate import ode
import pygame
import math

type TState = tuple[float, float, float, float]
type TController = Callable[[TState], float]
class TParameters(TypedDict):
    R: float
    L: float
    g: float
    mC: float
    mR: float
    mMW: float
    mW: float
    JR: float
    JW: float
    JMW: float
class TEnvironmentConfiguration(TypedDict):
    SCREEN_WIDTH: int
    SCREEN_HEIGHT: int
    FPS: int
    BG_COLOR: tuple[int, int, int]
    FLOOR_COLOR: tuple[int, int, int]
    CART_COLOR: tuple[int, int, int]
    FLOOR_Y: int
    CART_WIDTH_M: float
    CART_HEIGHT_M: float
    PX_TO_M: float
    DISTURBANCE_FORCE: float
class TEnvironmentConfigurationPartial(TypedDict, total=False):
    SCREEN_WIDTH: int
    SCREEN_HEIGHT: int
    FPS: int
    BG_COLOR: tuple[int, int, int]
    FLOOR_COLOR: tuple[int, int, int]
    CART_COLOR: tuple[int, int, int]
    FLOOR_Y: int
    CART_WIDTH_M: float
    CART_HEIGHT_M: float
    PX_TO_M: float
    DISTURBANCE_FORCE: float


_BASE_STATE: Final[TState] = (0.0, 0.0, 0.0, 0.0)
_BASE_CONFIG: Final[TEnvironmentConfiguration] = {
    'SCREEN_WIDTH': 960,
    'SCREEN_HEIGHT': 640,
    'FPS': 60,
    'BG_COLOR': (126, 126, 210),
    'FLOOR_COLOR': (0, 0, 0),
    'CART_COLOR': (255, 255, 255),
    'FLOOR_Y': 9 * 640 // 10,
    'CART_WIDTH_M': .5,
    'CART_HEIGHT_M': .25,
    'PX_TO_M': 400,
    'DISTURBANCE_FORCE': 2,
}


class Pendulum:
    """
    Pendulum on a motorized cart simulator.
    This class handles both the rendering and the simulation of the system.
    """

    _controller: TController | None = None
    _parameters: TParameters
    _state: TState
    _env_conf: TEnvironmentConfiguration
    _solver: ode
    _t: float

    _screen: pygame.Surface
    _clock: pygame.time.Clock

    @property
    def state(self) -> TState:
        """
        Get the current state of the system.

        @return: Current state as a tuple (x, dx, p, dp)
        """
        return self._state

    @property
    def parameters(self) -> TParameters:
        """
        Get the system parameters.

        @return: System parameters as a dictionary
        """
        return self._parameters
    
    @property
    def env_configuration(self) -> TEnvironmentConfiguration:
        """
        Get the environment configuration.

        @return: Environment configuration as a dictionary
        """
        return self._env_conf

    def _ddx(self, state: TState, T: float, U: float) -> float:
        x, dx, p, dp = state
        R = self._parameters['R']
        L = self._parameters['L']
        g = self._parameters['g']
        mC = self._parameters['mC']
        mR = self._parameters['mR']
        mMW = self._parameters['mMW']
        mW = self._parameters['mW']
        JR = self._parameters['JR']
        JW = self._parameters['JW']
        JMW = self._parameters['JMW']
        s1 = math.sin(p)
        s2 = math.sin(2*p)
        c2 = math.cos(2*p)
        return (R*(2*JR*U + 2*JR*R*T + 2*L**2*U*mR - 2*L**2*R*T*mR*c2 + 2*L**3*R*dp**2*mR**2*s1 - L**2*R*g*mR**2*s2 + 2*JR*L*R*dp**2*mR*s1))/(2*JMW*JR + 2*JR*JW + L**2*R**2*mR**2 + 2*JMW*L**2*mR + 2*JW*L**2*mR + 2*JR*R**2*mC + 2*JR*R**2*mMW + 2*JR*R**2*mR + 2*JR*R**2*mW + 2*L**2*R**2*mC*mR + 2*L**2*R**2*mMW*mR + 2*L**2*R**2*mR*mW - L**2*R**2*mR**2*c2)

    def _ddp(self, state: TState, T: float, U: float) -> float:
        x, dx, p, dp = state
        R = self._parameters['R']
        L = self._parameters['L']
        g = self._parameters['g']
        mC = self._parameters['mC']
        mR = self._parameters['mR']
        mMW = self._parameters['mMW']
        mW = self._parameters['mW']
        JR = self._parameters['JR']
        JW = self._parameters['JW']
        JMW = self._parameters['JMW']
        s1 = math.sin(p)
        s2 = math.sin(2*p)
        c1 = math.cos(p)
        return (L*(2*JMW*T*c1 + 2*JW*T*c1 + R**2*g*mR**2*s1 - R*U*mR*c1 + JMW*g*mR*s1 + JW*g*mR*s1 + 2*R**2*T*mC*c1 + 2*R**2*T*mMW*c1 + R**2*T*mR*c1 + 2*R**2*T*mW*c1 - (L*R**2*dp**2*mR**2*s2)/2 + R**2*g*mC*mR*s1 + R**2*g*mMW*mR*s1 + R**2*g*mR*mW*s1))/(JMW*JR + JR*JW + L**2*R**2*mR**2 + JMW*L**2*mR + JW*L**2*mR + JR*R**2*mC + JR*R**2*mMW + JR*R**2*mR + JR*R**2*mW + L**2*R**2*mC*mR + L**2*R**2*mMW*mR + L**2*R**2*mR*mW - L**2*R**2*mR**2*c1**2)

    def _make_rhs(self):
        """
        Computes the derivative of the state vector and returns it.
        """
        def rhs(t, state, T, U):
            return [state[1], self._ddx(state, T, U), state[3], self._ddp(state,T, U)]

        return rhs


    def __init__(
        self,
        parameters: TParameters,
        initial_state: TState | None = None,
        initial_time: float = 0.0,
        env_configuration: TEnvironmentConfigurationPartial | None = None,
    ) -> None:
        """
        Initialize the pendulum simulator.

        @param parameters: Dictionary containing the system parameters.
        Expected keys: 'R', 'L', 'g', 'mC', 'mR', 'mMW', 'mW', 'JR', 'JW', 'JMW'

        @param initial_state: List containing the initial state of the system.
        Expected format: [x, dx, p, dp]

        @param env_configuration: Dictionary containing the environment configuration.
        Expected keys: 'SCREEN_WIDTH', 'SCREEN_HEIGHT', 'FPS', 'BG_COLOR', 'FLOOR_COLOR', 'CART_COLOR', 'FLOOR_Y', 'CART_WIDTH_M', 'CART_HEIGHT_M', 'PX_TO_M', 'DISTURBANCE_FORCE'

        @return: None
        """
        self._parameters = parameters
        self._state = _BASE_STATE if initial_state is None else initial_state
        self._t = initial_time
        self._env_conf = _BASE_CONFIG.copy()
        if env_configuration is not None:
            self._env_conf.update(env_configuration)
        if self._t < 0.0:
            raise ValueError("Initial time must be non-negative.")
        self._solver = ode(self._make_rhs()).set_integrator('dopri5', rtol=1e-8, atol=1e-10)
        self._solver.set_initial_value(self._state, self._t)

        pygame.init()
        self._screen = pygame.display.set_mode((self._env_conf['SCREEN_WIDTH'], self._env_conf['SCREEN_HEIGHT']))
        pygame.display.set_caption("Pendulum on a Motorized Cart Simulation")
        self._clock = pygame.time.Clock()


    def _render(self) -> None:
        """
        Render the current state of the system on the screen.
        """
        self._screen.fill(self._env_conf['BG_COLOR'])
        pxm = self._env_conf['PX_TO_M']
        c = self._env_conf['SCREEN_WIDTH']//2 - self._state[0] * pxm
        wheel_r = self._parameters['R'] * pxm
        wheel_dia = 2 * wheel_r

        # Floor
        pygame.draw.rect(
            self._screen,
            self._env_conf['FLOOR_COLOR'],
            ((0, self._env_conf['FLOOR_Y']),
            (self._env_conf['SCREEN_WIDTH'], self._env_conf['SCREEN_HEIGHT'])),
            0
        )
        pygame.draw.line(
            self._screen,
            tuple(0xFF - c for c in self._env_conf['FLOOR_COLOR']),
            (self._env_conf['SCREEN_WIDTH']//2,self._env_conf['FLOOR_Y'] + 5),
            (self._env_conf['SCREEN_WIDTH']//2,self._env_conf['FLOOR_Y'] + 40),
            4
        )

        # Cart
        pygame.draw.rect(
            self._screen,
            self._env_conf['CART_COLOR'],
            ((c-self._env_conf['CART_WIDTH_M']*pxm/2, self._env_conf['FLOOR_Y'] - self._env_conf['CART_HEIGHT_M']*pxm - wheel_dia + 3),
            (self._env_conf['CART_WIDTH_M']*pxm,self._env_conf['CART_HEIGHT_M']*pxm - 5)),
            0
        )
        pygame.draw.line(
            self._screen,
            self._env_conf['CART_COLOR'],
            (c, self._env_conf['FLOOR_Y'] - min(20, wheel_dia*2//3)),
            (c, self._env_conf['FLOOR_Y'] - min(3, wheel_dia//3)),
            4
        )

        # Wheels
        c1 = (c-self._env_conf['CART_WIDTH_M']*pxm/2 + wheel_r, self._env_conf['FLOOR_Y'] - wheel_r)
        c2 = (c+self._env_conf['CART_WIDTH_M']*pxm/2 - wheel_r, self._env_conf['FLOOR_Y'] - wheel_r)
        pygame.draw.circle(
            self._screen,
            self._env_conf['CART_COLOR'],
            c1, wheel_r, 0
        )
        pygame.draw.line(
            self._screen,
            tuple(0xFF - c for c in self._env_conf['CART_COLOR']),
            c1, (c1[0] - 0.9*wheel_r * math.sin(self._state[0]/self._parameters['R']), c1[1] - 0.9*wheel_r * math.cos(self._state[0]/self._parameters['R'])),
            3
        )
        pygame.draw.circle(
            self._screen,
            self._env_conf['CART_COLOR'],
            c2, wheel_r, 0
        )
        pygame.draw.line(
            self._screen,
            tuple(0xFF - c for c in self._env_conf['CART_COLOR']),
            c2, (c2[0] - 0.9*wheel_r * math.sin(self._state[0]/self._parameters['R']), c2[1] - 0.9*wheel_r * math.cos(self._state[0]/self._parameters['R'])),
            3
        )

        # Rod
        pivot = (c, self._env_conf['FLOOR_Y'] - self._env_conf['CART_HEIGHT_M']*pxm - wheel_dia)
        pygame.draw.line(
            self._screen,
            tuple(0xFF - c for c in self._env_conf['CART_COLOR']),
            pivot,
            (pivot[0] - self._parameters['L'] * pxm * math.sin(self._state[2]), pivot[1] - self._parameters['L'] * pxm * math.cos(self._state[2])),
            6
        )
        pygame.draw.line(
            self._screen,
            self._env_conf['CART_COLOR'],
            (pivot[0] - 2*math.sin(self._state[2]), pivot[1] - 2*math.cos(self._state[2])),
            (pivot[0] - (self._parameters['L'] * pxm -2) * math.sin(self._state[2]), pivot[1] - (self._parameters['L'] * pxm -2) * math.cos(self._state[2])),
            3
        )

        pygame.display.flip()


    def _step(self, dt: float, T: float, U: float) -> None:
        """
        Performs a single step of the simulation.
        Time elapsed from the last step is dt.
        The control inputs T and U are assumed constant during this step.
        """
        self._solver.set_f_params(T, U)
        self._solver.integrate(self._t + dt)
        if not self._solver.successful():
            raise RuntimeError("Integration failed.")
        self._t = self._solver.t
        self._state = tuple(self._solver.y.copy())
    

    def set_controller(self, controller: TController) -> None:
        """
        Set the controller for the system.

        @param controller: A callable that takes the current state and returns a control input.
        """
        self._controller = controller


    def run(self) -> None:
        """
        Run the simulation and rendering loop.
        """
        fps = self._env_conf['FPS']
        iT = 0.0
        iU = 0.0

        while True:
            dt = self._clock.tick(fps) / 1000.0

            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    pygame.quit()
                    return
                if event.type == pygame.KEYDOWN and event.key == pygame.K_ESCAPE:
                    pygame.quit()
                    return
                if event.type == pygame.KEYDOWN and event.key == pygame.K_a:
                    iU = self._env_conf['DISTURBANCE_FORCE']
                if event.type == pygame.KEYDOWN and event.key == pygame.K_d:
                    iU = -self._env_conf['DISTURBANCE_FORCE']
                if event.type == pygame.KEYUP and event.key in [pygame.K_a, pygame.K_d]:
                    iU = 0.0
            
            if self._controller is not None:
                iT = self._controller(self._state)

            self._step(dt, iU, iT)
            self._render()