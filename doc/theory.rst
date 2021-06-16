******
Theory
******
Solar
=====
Solar Angles
------------
While the output power of the Sun is usually considered as a constant, the
amount of power arriving at the Earth's surface varies according to the time,
location, weather, and relative position of the Earth with respect to the Sun.
Besides, the available data needed for a solar power calculation is usually given
for a horizontal surface and most of the PV systems are placed in a tilted
position. Therefore, it is necessary to calculate a set of parameters describing
the Sun's relative position with respect to the position of the system being
irradiated.
These parameters are calculated for points located at the center of every
pixel (with high resolution) within the extension under analysis and for every
hour of the year.

Declination Angle :math:`\delta`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This angle varies during the year due to the tilt of the 
Earth's axis, which is 23.45° tilted, so the declination ranges between -23.4°
and 23.45° through the year. The declination could be interpreted as the latitude
where the Sun's rays perpendicularly strike the Earth's surface at solar noon.
This value is the same for all the locations within the globe for a given day and
is calculated as follows :cite:`Scharmer.2000`:

.. math::

	\delta=\arcsin\left(0.3978 \sin\Big(\frac{2\pi N}{365.25}-1.4+0.0355 \sin\Big(\frac{2\pi N}{365.25}-0.0489\Big)\Big)\right)

where N is the day of the year.

Solar Time
^^^^^^^^^^
The time is important to define the position of the Sun in the
sky. However, it is easier to use the time if it is converted into solar time. To
do so, a few corrections are needed. The equation of time is an empirical
equation which corrects the error caused by the axial tilt of the Earth and the
eccentricity of its orbit :cite:`Scharmer.2000`:

.. math::
	EOT=-0.128 \sin\Big(\frac{360}{365.25}N-2.8\Big)-0.165 \sin\Big(\frac{720}{365.25}N +19.7\Big)

When the time is given in GMT, as it is for this model, it is also necessary to
take into account the longitude of the location, hence the time correction:

.. math::
	TC=EOT+longitude/15

where the factor 15 accounts for the geographical span of each time zone (15°
of longitude). With this correction, the local solar time is calculated as:

.. math::
	LST = T_{GMT}+TC

where LST is the local solar time. An even more apporpriate time measure 
for solar calculations is the hour angle :math:`\omega` This converts hours into 
degrees which indicate how the Sun moves in the sky relatively to the Earth, 
where the solar noon is 0°, the anlges after the noon positive, and before 
the noon negative.

.. math::
	\omega=15 (LST-12)

Another important quantity is the duration of the day, which is delimited by the
sunrise and the sunset. The sunrise and sunset have the same value, however,
the sunrise is considered negative and the sunset positive. They depend on the
day of the year and the location on the Earth (denoted by the declination and
the latitude :math:`\phi` respectively) and they are calculated for a horizontal surface as
follows :cite:`Klein.1977`:

.. math::
	\omega_s=\arccos(-\tan\phi\tan\delta)
	
Due to self-shading, a tilted plane might be exposed to different sunrise and
sunset's values. Also, if the plane is not facing the equator, the sunrise and
sunset angles will be numerically different for such surfaces. The following
equations consider this orientation changes for the sunrise and sunset values of
a tilted plane :cite:`Klein.1977`.

.. math::
	a=\frac{\cos\phi}{\tan\beta}+\sin\phi
.. math::
	b=\tan\delta\cos\phi\cos\gamma-\frac{\sin\phi}{\tan\beta}
.. math::
	\omega_s'=\cos\Big[\frac{ab\pm\sin\gamma\sqrt{a^2-b^2+\sin^2\gamma}}{a^2+\sin^2\gamma}\Big]

where :math:`\gamma` is the azimuthal orientation of the panel and :math:`\beta` is the tilt of the panel
(for this model, chosen as the optimal tilt according to the latitude). These
equations might give higher values than the real sunrise and sunset values. This
would imply that the Sun rises over the tilted plane before it has risen over the
horizon or that when the Sun sets, there is still light striking the plane. As this
is wrong, the sunrise and sunset values for a horizontal plane must be compared
with the values for a tilted plane and the lower values (for both sunrise and
sunset) must be selected.

.. math::
	\omega_0=\min(\omega_s,\omega_s')

Incidence Angle :math:`\theta`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As stated before, PV panels are not normally parallel to
the Earth's surface, so it is necessary to calculate the incidence angle of the
Sun's rays striking the surface of the panel. Nevertheless, a set of angles must
be calculated in order to calculate the incidence angle.

The elevation angle :math:`\alpha` or altitude angle measures the angular distance
between the Sun and the horizon. It ranges from 0° at the sunrise to 90° at the
noon (the value at the noon varies depending on the day of the year) :cite:`William.2014`.

.. math::
	\alpha=\arcsin[\sin\delta \sin\phi+\cos\delta\cos\phi\cos\omega]

The azimuth angle :math:`A_{z}` is an angular measurement of the horizontal position
of the Sun. It could be seen as a compass direction with 0° to the North and
180° to the South. The range of values of the azimuth angle varies over the
year, going from 90° at the sunrise to 270° at the sunset during the equinoxes.
The equation for the azimuth depends on the time of the day. For the solar
morning, it is :cite:`William.2014`:

.. math::
	Az_{am}=\arccos\Big(\frac{\sin\delta\cos\phi-\cos\delta\sin\phi\cos\omega}{\cos\alpha}\Big)

and for the afternoon:

.. math::
	Az_{pm}=360-Az_{am}

With the already calculated angles, it is possible to calculate the incidence angle,
which is the angle between the surface's normal and the Sun's beam radiation :cite:`Reindl.1988`:

.. math::
	\begin{split}
	\theta_{i}=&\arccos(\sin\delta \sin\phi \cos\beta \\
	&-\sin\delta\cos\phi\sin\beta\cos\gamma \\
	&+\cos\delta\cos\phi\cos\beta\cos\omega \\
	&+\cos\delta\sin\phi\sin\beta\cos\gamma\cos\omega\\
	&+\cos\delta\sin\beta\sin\gamma\sin\omega)
	\end{split}
	
Tracking
^^^^^^^^
When one-axis tracking is active, the tilt angle :math:`\beta` and the azimuthal orientation :math:`\gamma` 
of the panel change constantly as the panel follows the sun. In this model a tilted on-axis tracking with east-west tracking is considered.
The rotation of the plane around the axis is deffned by the rotation angle
R, it is calculated in order to achieve the smallest incidence angle for the plane
by the following equations :cite:`William.2013`:

.. math::
	X=\frac{-\cos\alpha \sin(A_z-\gamma_a)}{-\cos\alpha \cos(A_z-\gamma_a)\sin\beta_a+\sin\alpha \cos\beta_a}
.. math::
	\Psi=
    \begin{cases}
      0, & \text{if}\ X=0, \text{ or if } X>0 \land (A_z-\gamma_a)>0, \text{ or if } X<0 \land (A_z-\gamma_a)<0 \\
      180, & \text{if}\ X<0 \land (A_z-\gamma_a)>0 \\
      -180, & \text{if}\ X>0 \land (A_z-\gamma_a)<0
    \end{cases}


for the previous equations, :math:`\beta_a` and :math:`\gamma_a` are considered as the tilt and azimuthal
orientation of the tracking axis respectively. The variable :math:`\Psi` places R in the correct 
trigonometric quadrant. For the selection of :math:`\Psi`, the difference :math:`(A_z-\gamma_a)` must be considered 
as the angular displacement with the result within the range
of -180°to 180°. Once the rotation angle is calculated, the tilt and azimuthal
orientation of the panel are calculated as follows:

.. math::
	\beta=\arccos(\cos R \cos\beta_a)
.. math::
	\gamma=
    \begin{cases}
      \gamma_a+\arcsin\Big(\dfrac{\sin R}{\sin\beta}\Big), & \text{for}\ \beta_a \neq 0, -90 \leq R \leq 90 \\
      \gamma_a-180-\arcsin\Big(\dfrac{\sin R}{\sin\beta}\Big), & \text{for}\ -180 \leq R < -90 \\
      \gamma_a+180-\arcsin\Big(\dfrac{\sin R}{\sin\beta}\Big), & \text{for}\ 90 < R \leq -90 
    \end{cases}

Then the incidence angle is calculated using the new :math:`\beta` and :math:`\gamma` angles. For two-axis tracking the beta angle is considered as the complementary
angle to the altitude angle while the azimuthal orientation angle :math:`\gamma` is considered as :math:`A_z - 180`.

Solar Power
-----------
To calculate the solar power we must start with the solar constant. However,
the irradiance striking the top of the Earth's atmosphere (TOA) varies over the
year. This is due to the eccentricity of the Earth's orbit and its tilted axis. The
TOA is calculated according to the following equation :cite:`Reindl.1988`:

.. math::
	TOA=G_{sc}\Big[1+0.03344  \cos\Big(\frac{2 \pi N}{365.25}-0.048869\Big)\Big] \sin\alpha

where :math:`G_{sc}` is the solar constant (1367 W/m2) and N is the day of the year.
This equation escalates the solar constant by multiplying it with a factor related
to the eccentricity of the Earth's orbit and the sinus of the altitude angle of the
Sun to finally get the extraterrestrial horizontal radiation.

To calculate the amount of horizontal radiation at the surface of the Earth,
the attenuation caused by the atmosphere must be considered. A way to measure
this attenuation is by using the clearness index :math:`k_t`. This value is the ratio
between the extraterrestrial horizontal radiation and the radiation striking the
Earth's surface:

.. math::
	k_t=\frac{GHI_{M2}}{TOA_{M2}}

where :math:`GHI_{M2}` and :math:`TOA_{M2}` are the global horizontal irradiance and the top
of the atmosphere radiation extracted from MERRA-2 data. Furthermore, the
GHI is made-up by diffuse and beam radiation:

.. math::
	GHI=G_b+G_d

where :math:`G_b` is the beam radiation, which is the solar radiation that travels directly
to the Earth's surface without any scattering in the atmosphere, and :math:`G_d` stands
forfor the diffuse radiation, the radiation that comes to a surface from all directions
as its trajectory is changed by the atmosphere. These two components have
different contributions to the total irradiance on a tilted surface, so it is necessary
to distinguish between them. This can be done using the correlation of Erbs
et al :cite:`Erbs.1982`, which calculates the ratio R of the beam and diffuse radiation as a
function of the clearness index.

.. math::
	R=
    \begin{cases}
      1 - 0.09 k_t, & \text{for}\ k_t\leq0.22 \\
      0.9511 - 0.1604 k_t+ 4.388 k_t^2 - 16.638 k_t^3 + 12.336 k_t^4, & \text{for}\ 0.22> k_t \leq 0.8 \\
      0.165, & \text{for}\ k_t>0.8
    \end{cases}
	
Furthermore, diffuse radiation could be divided into more components. The
HDKR model :cite:`Duffie.2013,Reindl.1988`,  developed by Hay, Davies, Klucher, and Reindl in 1979
assumes isotropic diffuse radiation, which means that the diffuse radiation is
uniformly distributed across the sky. However, it also considers a higher radiation
intensity around the Sun, the circumsolar diffuse radiation, and a horizontal
brightening correction. To use the HDKR model, some factors must be defined
first :cite:`Reindl.1988`. The ratio of incident beam to horizontal beam:

.. math::
	R_b=\frac{\cos\theta_i}{\sin\alpha}
	
The anisotropy index for forward scattering circumsolar diffuse irradiance:

.. math::
	A_i=(1-R)k_t
	
The modulating factor for horizontal brightening correction:

.. math::
	f=\sqrt{1-R}

Then the total radiation incident on the surface is calculated with the next equations:

.. math::
	GHI=k_tTOA
.. math::
	\begin{split}
    G_T=GHI\Big[(1-R+RA_i)R_b+R(1-A_i)\Big(\frac{1+\cos\beta}{2}\Big)&\Big(1+f\sin^3\frac{\beta}{2}\Big)\\
    &+\rho_g\Big(\frac{1-\cos\beta}{2}\Big)\Big]
	\end{split}

Where :math:`\rho_g` is the ground reflectance or albedo and it is related to 
the land use type of the location under analysis. The first term of this equation 
corresponds to the beam and circumsolar diffuse radiation, the second to the isotropic 
and horizon brightening radiation, and the last one to the incident ground-reflected radiation. 

PV
--
Temperature losses
^^^^^^^^^^^^^^^^^^
In a PV panel, not all the radiation absorbed is converted
into current. Some of this radiation is dissipated into heat. Solar cells, like all
other semiconductors, are sensitive to temperature. An increase of temperature
results in a reduction of the band gap of the solar cell which is translated into
a reduction of the open circuit voltage. The overall effect is a reduction of the
power output of the PV system. To calculate the power loss of a solar cell it is
necessary to know its temperature. This can be expressed as a function of the
incident radiation and the ambient temperature :cite:`Maturi.2014`:

.. math::
	T_{cell}=T_{amb}+kG_T

where :math:`T_{amb}` is the ambient temperature and k is the Ross coefficient, which
depends on the characteristics related to the module and its environment. It
is defined based on the land use type of the region where the panel is located.
With the temperature of the panel, the fraction of the irradiated power which
is lost can be calculated as:

.. math::
	Loss_T=(T_{cell}-T_r)T_k

where :math:`T_r` is the rated temperature of the module according to standard test
conditions and :math:`T_k` is the heat loss coefficient. Both values are usually given on
the data sheets of the PV modules.

PV Capacity Factor Calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
It is the ratio of the actual power output to the theoretical
maximum output which is normally considered as :math:`1000 W/m^{2}`. The temperature
loss is also considered for this calculation:

.. math::
	CF_{PV} = \frac{G_T(1-Loss_T)}{1000}
	
Ground coverage ratio (GCR)
^^^^^^^^^^^^^^^^^^^^^^^^^^^
It is also important to consider the area lost due to the space between the modules 
or due to the modules shading adjacent modules. This is done with the GCR which is 
the ratio of the module area to the total ground area.

.. math::
	GCR=\frac{1}{\cos\beta+|\cos A_z| \cdot \Big(\dfrac{\sin\beta}{\tan \alpha}\Big)}

CSP
---
For its popularity and long development history, the parabolic trough technology was chosen to model Concentrated Solar power.

Convection Losses
^^^^^^^^^^^^^^^^^
The receiver of parabolic troughs are kept in a vacuum glass tube to prevent convection as much as possible. 
Radiative heat losses are still present and ultimatly results in convective losses between the glass tube and the air.
These heat losses are increased when wind is blowing around the receiver. The typical heat losses for a receiver 
can be estimated through the following empirical equation :cite:`Duffie.2013`:

.. math::
	Q_{Loss} =  A_r(U_{L_{cst}} + U_{L_{Wind}} \cdot V_{Wind}^{0.6})(T_i-T_a)

where :math:`A_r` is the outer area of the receiver, :math:`U_{L_{cst}}` correspond to a loss coefficient at zero wind speed, 
:math:`U_{L_{Wind}}` is a loss coefficient dependent on the wind speed :math:`V_{Wind}`, 
:math:`T_i` is the average heat transfer fluid temperature, and :math:`T_a` is the ambient temperature.

Typical values for the :math:`U_{L_{cst}}` and :math:`U_{L_{Wind}}` are 1.06 :math:`kW/m^{2}K` and 1.19 :math:`kW/(m^{2}K(m/s)^{0.6})` respectively

Flow Losses
^^^^^^^^^^^
Flow loss coefficient or heat removal factor :math:`F_r` is the ratio between the actual heat transfer to the maximum heat transfer possible between 
the receiver and the heat transfer fluid (HTF). These losses result from the difference between the temperature of the receiver and 
the temperature of the HTF and are dependent on the heat capacity and the flow rate of the HTF. 
A typical value for parabolic troughs is 95%. 

CSP Capacity Factor Calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The capacity factor of a solar field is the ratio of the actual useful heat collected to the theoretical maximum heat output of 1000 W/m². It is given by the formula:

.. math::
	CF_{csp} = \frac{F_r(S - Q_{Loss})}{1000}
	
Where :math:`S` is the component of the DNI captured by the collector at an angle (based on one axis traking), 
:math:`Q_{Loss}` is the heat convection losses, and :math:`F_r` is the heat removal factor.

Wind
====
Wind Speed
----------
In order to use a single value for the following wind power calculations, a norm is calculated 
as if both variables were vectors, so the absolute velocity is:

.. math::
	W_{50M} = \sqrt{V50M^2+U50M^2}

However, the wind velocities extracted from MERRA-2 are given for a height
of 50 meters, which does not correspond to the hub height of the wind turbines;
therefore, they must be extrapolated.

Wind Shear
----------
While the wind is hardly affected by the Earth's surface at a
height of about one kilometer, at lower heights in the atmosphere the friction of
the Earth's surface reduces the speed of the wind :cite:`Danish.2003`. One of the most common
expressions describing this phenomenon is the Hellmann exponential law, which
correlates the wind speed at two different heights :cite:`Kaltschmitt.2007`.

.. math::
	v=v_0\Big(\frac{H}{H_0}\Big)^\alpha

Where :math:`v` is the wind speed at a height :math:`H`, :math:`v_0` is the wind speed at a height :math:`H_0`
and :math:`\alpha` is the Hellmann coefficient, which is a function of the topography and air
stability at a specific location.

Wind Power
----------
The wind turbines convert the kinetic energy of the air into torque. The
power that a turbine can extract from the wind is described by the following
expression :cite:`Masters.2004`:

.. math::
	P=\frac{1}{2}\rho A v^3 C_p

where :math:`\rho` is the density of the air, :math:`v` is the speed of the wind and :math:`C_p` is the power
cofficient. As it is shown in the previous equation, the energy in the wind varies
proportionally to the cube of the wind's speed. Therefore, the power output of
wind turbines is normally described with cubic power curves. However, there
are some regions within those curves which have special considerations.

Cut-in wind speed
^^^^^^^^^^^^^^^^^
The wind turbines start running at a wind speed between
3 and 5 m/s to promote torque and acceleration. Therefore, there is no power
generation before this velocity.

Rated Wind Speed
^^^^^^^^^^^^^^^^
The power output of a wind turbine rises with the wind
speed until the power output reaches a limit defined by the characteristics of
the electric generator. Beyond this wind speed, the design and the controllers of
the wind turbine limit the power output so this does not increase with further
increases of wind speed.

Cut-out wind speed
^^^^^^^^^^^^^^^^^^
The wind turbines are programmed to stop at wind
velocities above their rated maximum wind speed(usually around 25 m/s) 
to avoid damage to the turbine and its surroundings,
so there is no power generation after this velocity.

Wind Onshore and Offshore Capacity factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Finally, the capacity factors, which are the ratios of the
actual power output to the theoretical maximum output (rated wind speed),
are calculated according to the previously presented regions:

.. math::
	CF=
    \begin{cases}
      \dfrac{W_{hub}^3-W_{in}^3}{W_r^3-W_{in}^3}, & \text{for}\ W_{in}<W_{hub}<W_r \\
      1, & \text{for}\ W_r\leq W_{hub}\leq W_{out} \\
      0, & \text{for}\ W_{hub}<W_{in} | W_{hub}>W_{out}
    \end{cases}

where :math:`W_{in}`Win is the cut-in wind speed, :math:`W_{out}` is the cut-out wind speed, and :math:`W_r`
is the rated wind speed. The area is not included in the previous equations as
it does not change in both generation states (actual and theoretical maximum
power). While the density could vary for both states, the overall impact of a
change in density is negligible compared to the wind speed and therefore is
not included in the calculation.
