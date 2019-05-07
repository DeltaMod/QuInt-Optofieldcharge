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
+ Figures: Contains Misc. Images from earlier versions of the code 
------------------------------------------------------------------------------------------------------------------
+ Functions: Contains old versions of fun_Q and fun_QDelt that calculate the optical-field-induced currents
------------------------------------------------------------------------------------------------------------------
+ SavedPresets: Contains presets that can be loaded into the GUI to run a specific plot
----- + FigurePlotting: Presets used for the thesis. THESISGRAPH is a variable found inside UIGraphPlotting, it is 
                        used for more complicated graphs - e.g. Figure 14 from the thesis
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
