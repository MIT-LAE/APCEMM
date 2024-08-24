This folder contains python scripts to produce idealised met files containing trapezoidal moist layers (in an altitude
vs RHi plot).

The trapezoid is symmetrical about its centre, and is parametrised as follows:
	- centre_altitude_m: The centre altitude in meters.
	- moist_depth_m: The moist layer depth in meters; the vertical extent of the moist layer.
	- RHi_background_PC: The background RHi in %, which is also the lowest RHi value in the trapezium.
	- RHi_peak_PC: The peak trapezium RHi in %. If too shallow a gradient is used, this value might not be reached.
	- grad_RHi_PC_per_m: The gradient of the trapezium sides in RHi % / vertical m. If too shallow a gradient is used, 
							the moist layer will be a triangle instead of a trapezium and its peak will be lower than 
							or equal to the RHi_peak_PC value.
							
The met generator has additional parameters:
	- alt_resolution_m: The vertical interval (in m) at which the met is sampled at.
	- upper_alt_m: The vertical upper limit (in m) for met sampling.
	- shear: The horizontal wind shear (in 1/s). Constant throughout the entire met domain.
	- T_offset_K: A temperature offset (in K) from the ISA temperature distribution. Constant throughout the entire met
					domain.

Other than the RHi and shear, all other weather parameters assume the use of the International Standard Atmosphere.
The temperature is decoupled from the pressure by assuming a constant alpha regardless of any temperature offset.
The ideal gas relation is then used to calculate the air density from pressure and temperature to guarantee
consistency.

The python version used to develop this code is python 3.9.15, with the following modules:
	- matplotlib 3.6.2
	- netcdf4 1.6.3
	- numpy 1.23.5
	- xarray 2023.12.0
	
Desired improvements:
	- No distinction is made between dry and moist air density in this code. The difference is very small at altitude,
		but it is still a deficiency that needs to be addressed.
		
Contact:
	- Please contact ca525@cam.ac.uk (or Calebsakhtar on GitHub) for any queries or suggestions.
	
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.