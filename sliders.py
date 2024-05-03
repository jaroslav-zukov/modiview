import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

class DiscreteSlider(Slider):
    """A matplotlib slider widget with discrete steps."""
    def set_val(self, val):
        discrete_val = round(val)
        super(DiscreteSlider, self).set_val(discrete_val)

# Initial parameters
amp = 5
freq = 3

# Create the figure and the line that we will manipulate
fig, ax = plt.subplots()
line, = plt.plot(np.linspace(0, 1, 1000), np.sin(2 * np.pi * freq * np.linspace(0, 1, 1000)))

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

# Create the slider
sfreq = DiscreteSlider(axfreq, 'Freq', 0.1, 30.0, valinit=freq)

# Update function that is called whenever the slider value is changed
def update(val):
    freq = sfreq.val  # The slider value is already an integer
    line.set_ydata(np.sin(2 * np.pi * freq * np.linspace(0, 1, 1000)))
    fig.canvas.draw_idle()

# Call the update function whenever the slider value is changed
sfreq.on_changed(update)

plt.show()