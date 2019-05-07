# QuInt-Optofieldcharge
Reproduction and expansion of a models created by Jacob B. Khurgin. Includes a GUI that allows for the plotting of variable parameter laser pulses, and have the photoinduced charge on a nanodevice estimated from their dispersed components.
<pre>                                                                       
              __                                                    __ 
           __/ /___________________________________________________/ /_
          /_  __/_____/_____/_____/_____/_____/_____/_____/_____/_  __/
           ///    ______            __             __            ///   
          / /    / ____/___  ____  / /____  ____  / /______     / /    
         / /    / /   / __ \/ __ \/ __/ _ \/ __ \/ __/ ___/    / /     
        / /    / /___/ /_/ / / / / /_/  __/ / / / /_(__  )    / /      
       / /     \____/\____/_/ /_/\__/\___/_/ /_/\__/____/    / /       
      /_/ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ /_/        
     ( |_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_| )         
 _ _ |/_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ |/_ _ _     
( | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | )    
|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/|/   
 </pre>
 
 + KhurginUI.m: To launch the UI, run this script - It will open a window containing a figure plotter and variable selector.
			-Select a simulation with the dropdown, then select a plot. If you have run the simulation once, you can simply change the plots without needing to run the simulation again.
			-Click on Popout/Dock (The text is buggy, but the feature is not) graph to plot inside of a figure window - This will let you save your plots
			-Save Preset will save all currently entered variables (as well as the corresponding simulation and graph selection)
			-Load Preset will load the saved variables from a file
			-Everything else should be self explanatory :)
----------------------------------------------------------------------------------------------------------------------------------------+KhurginUI.fig - Needlessly oversized guide file which defines the UI - You just need it in the folder when running KhurginUI, and it should work. To make changes, type `guide` after setting your folder destination to this folder.
----------------------------------------------------------------------------------------------------------------------------------------
+ FigPos.m - Small script that detects your screensize and divides it into a 3x3 figure plotting grid. 
Use: f = figure(1); f.Position = FigPos(1,1). Documentation exists in the script, `help FigPos`
----------------------------------------------------------------------------------------------------------------------------------------
+ PhotoinducedChargetoUI.m - Script responsible for runnig a fourier limited version of the code. It calls fun_Q after determining a(t).
			-Note, most of its functionality can be replicated (and better) in DispersionToUI.m
---------------------------------------------------------------------------------------------------------------------------------------- 
+ DispersionToUI.m - A script that handles executing the functions necessary to calculate the temporal pulse, to fourier transform said pulse and apply dispersion, before transforming it back. It also handles the temporal delay implementation for the two colour, and cross polarised case It does not require modification, and while it is possible to run as standalone, is designed to be run from the UI
		- NOTE: Make sure to only set n_min/n_max greater than [-1,1] if you do intend to plot a delay - otherwise, you will be waiting for a long time
--------------------------------------------------------------------------------------------------------------------------------------- 
+ Variables.m - Sets standard variables and calculates things like the refractive index (only from SiO2 right now), then takes dn/dlambda (first three terms), to find the terms that give determine the GD, GVD and TOD
----------------------------------------------------------------------------------------------------------------------------------------
+FFTD.m: A function that handles the fourier transform, dispersion application, and inverse fourier transform. Extensive documentation exists in the script itself, either that or call `help FFTD` to find out exactly how it works. 
---------------------------------------------------------------------------------------------------------------------------------------- 
+ Figures: Contains Misc. Images from earlier versions of the code 
------------------------------------------------------------------------------------------------------------------
+ Functions: Contains old versions of fun_Q and fun_QDelt that calculate the optical-field-induced currents
------------------------------------------------------------------------------------------------------------------
+ SavedPresets: Contains presets that can be loaded into the GUI to run a specific plot
    - FigurePlotting: Presets used for the thesis. THESISGRAPH is a variable found inside UIGraphPlotting, it is 
                        used for more complicated graphs - e.g. Figure 14 from the thesis. To reproduce, some results, 
                        you will need to find them (searching for THESIS will do it)
------------------------------------------------------------------------------------------------------------------
+ TinyFunctions - Contains fun_Q and fun_QDelt - they are used to calculate the photoinduced charge taking into 
                   consideration N number of terms from <a^2n+1>. Increasing ORD will increase this number.
==================================================================================================================
+ fun_Q(F_0x,F_a,a2disp,Aeff,Order): 

                                     F_0x   - Field strength (A*10^10 Vm^-1)      [Vector]
                                     F_a    - Atomic field strength (5.36*10^10)  [Single Term]
                                     a2disp - Integral under a(t)^2n+1            [Vector]
                                     Aeff   - Effective Area (10^-12 m^2)         [Single Term]
                                     Order  - Terms 2n+1 taken into consideration [Single Term](must have same n# terms of a2disp)
 
 Usage: Input terms (must have eps_0 =  8.85418782 * 10^-12 defined in body) get data.
 Output:  avecFN - stores string data as a2n - used to name ATermsQ, use to get data from ATermsQ.(avecFN{n})
          atFN   - stores string data of each term a(t)^2n_1   - useful for legends
          ATermsQ.(avecFN{n}) - stores the contributions to charge for each term
 Q = fun_Q 
+ fun_QDelt(F_0x,F_0y,F_a,a2disp,Aeff,Order)
                                      F_0y  - Injection field (assume F_0x is driving)
                                      

DispersionToUI.m	Add files via upload	11 hours ago
FFTD.m	Add files via upload	11 hours ago
FigPos.m	Add files via upload	18 days ago
KhurginUI.m	Add files via upload	8 days ago
PhotoInducedChargetoUI.m	Add files via upload	a day ago
PowerDependence.m	Add files via upload	4 days ago
ProcessPlotting.psd	Add files via upload	4 days ago
UIConsoleOutput.m	Add files via upload	18 days ago
UIGraphPlotter.asv	Add files via upload	4 days ago
UIGraphPlotter.m	Add files via upload	11 hours ago
UINCleaner.m	Add files via upload	18 days ago
UIReinitialise.m	Add files via upload	4 days ago
Untitled.m	Add files via upload	18 days ago
Variables.m
 
 Will change later
