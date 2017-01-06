===============            ====================            ==================
README.TXT File            SAP2000 Version 7.40            Release 2000/06/12
===============            ====================            ==================


PLEASE READ THIS FILE!

It contains important information that is more current than the Manuals.



=================
TABLE OF CONTENTS
=================

1. SAP2000 Version 7.40
2. CHANGES TO INSTALLATION INSTRUCTIONS
3. New Hardware Key Device for Nonlinear Version 7
4. Tutorial files on the Web
5. Significant Changes from Previous Versions
6. Upgrading to SAP2000 from SAP90 Version 5.4 or 5.5
7. Analysis Executables




=======================
1. SAP2000 Version 7.40
=======================

SAP2000 is a fully graphical Windows program, and is available in three
different versions: SAP2000 Standard, SAP2000 PLUS, and SAP2000 Nonlinear.
All versions require Windows 95/98/NT/2000 to be running on your computer.
Model creation, analysis, display of results, and design are all performed
within the SAP2000 graphical user interface.

SAP2000 is a very extensive program, featuring a friendly and powerful
graphical user interface and sophisticated analytical capabilities. However, a
few of the advanced analysis features cannot be edited in the graphical user
interface, such as the Plane, Asolid, and Solid elements. The geometry and
analysis results for these elements CAN be viewed in the graphical user
interface.

Users who have need of advanced analysis features described in the Analysis
Reference manual that cannot be edited in the graphical user interface can use
the input data text file as described in the Input Data File manual. See Topics
"Overview" and "Input Data Files and the Graphical User Interface" in Chapter
2 of the Input Data File manual for more information.

For users who already have SAP2000, we recommend that you remove previous
versions of the program and begin using version 7.40 as soon as possible in
order to take advantage of the latest features, improvements, and error
corrections.




=======================================
2. CHANGES TO INSTALLATION INSTRUCTIONS
=======================================

For a new installation from the CD, please read topic "Installation" in chapter
"Getting Started" in the "SAP2000 Getting Started" manual, and note the
following changes:

- To perform a "Single User" installation, run

     \Setup.exe

  on the SAP2000 CD. This installs the entire SAP2000 program on a local
  machine, ready to be used.

- To perform a "Network Server" installation, run

     \Network Server\Setup.exe

  on the SAP2000 CD. This installs the Server version of SAP2000 on a network
  server. The program cannot be run yet. Each workstation that will use SAP2000
  must perform a "Network Workstation" installation described next.

- To perform a "Network Workstation" installation, run

     \Setup.exe

  from the folder on the server where the "Network Server" installation of
  SAP2000 was performed, as described above. This installs a limited number of
  files on the local workstation, and makes SAP2000 ready to be used. The
  workstation must have access to the server installation whenever SAP2000 is
  used.

- The default installation folders for the Standard version are

     c:\Computers and Structures\SAP2000
     c:\Computers and Structures\SAP2000 Server
     c:\Computers and Structures\SAP2000 Workstation

- The default installation folders for the PLUS version are

     c:\Computers and Structures\SAP2000 Plus
     c:\Computers and Structures\SAP2000 Plus Server
     c:\Computers and Structures\SAP2000 Plus Workstation

- The default installation folders for the Nonlinear version are

     c:\Computers and Structures\SAP2000 Nonlinear
     c:\Computers and Structures\SAP2000 Nonlinear Server
     c:\Computers and Structures\SAP2000 Nonlinear Workstation

- For each of the three versions above:
     The first folder name is used for a "Single User" installation,
     The second folder name is used for a "Network Server" installation, and
     The third folder name is used for a "Network Workstation" installation.

- During a "Single User" or "Network Server" installation you will be asked
  for the program serial number. This can be found on the SAP2000 CD. This
  serial number is automatically entered when performing a "Network
  Workstation" installation. When running SAP2000, you can view the serial
  number that you entered from the Help > About SAP2000 menu item. This
  serial number may be needed for obtaining technical support and program
  upgrades.

- At the end of the installation process you will be given the option to have
  the Sentinel Driver installed. You should choose this option whenever you
  are installing on a machine that may be used to hold the hardware key device.
  You may also manually install the Sentinel Driver at any time by running

     KeyDrive\Setup.exe

  from the SAP2000 CD or from the SAP2000 folder after installation.

- When installing on a Windows 2000 machine, make sure that you have
  administrative rights.

- Running any SAP2000 Setup.exe a second time uninstalls the program. You
  will be given an opportunity to Cancel if you do this accidentally.




================================================
3. New Hardware Key Device for Nonlinear Version
================================================

SAP2000 Nonlinear version 7.xx uses a hardware key device that is different
from previous versions. If you are upgrading SAP2000 Nonlinear from a
version prior to 7.00, please return the old hardware key device to us as
soon as you receive this package.

The Standard and PLUS version of SAP2000 do not require a new hardware key
device.




============================
4. Tutorial files on the Web
============================

Tutorial files are available on our Web site in Adobe Acrobat Reader format.
You may download these files and follow the tutorials to learn more about the
use of SAP2000 and the new Static Pushover Analysis features of the program.

Please check our Web site occasionally for the presence of newer files. They
may be updated from time to time. Our Web address is:  www.csiberkeley.com




=============================================
5. Significant Changes from Previous Versions
=============================================

Listed in order, starting with the latest version:


Significant Changes from Version 7.21 to 7.40
---------------------------------------------

Improvements have been made to the "MS Access Database export" capability that
was added in version 7.21. This feature is still an early release and is being
tested and documented before final release.

The Ritz-vector algorithm has been modified so that purely static (massless)
modes can be calculated. These will be reported in the .OUT file as having a
period of 1.0E-101. Thus it is no longer strictly required that mass be
present at all degrees of freedom that are loaded by the Ritz starting load
vectors. However, including mass at loaded degrees of freedom is still
recommended whenever practical, and the Ritz-vector algorithm will still warn
you if it detects loads acting on massless degrees of freedom. These massless
modes are used in time-history analysis but have no significant effect in
response-spectrum analysis. Also, the orthogonalization procedure in the
Ritz-vector algorithm has been improved.

P-delta analysis now includes P-delta effects in all element types, not just
in Frame elements. For all elements except the Frame, P-delta effects are
based on the rigid-body rotations of the element and do not consider the
internal deformation of the element. The Nllink element balances the P-delta
moment with a shear couple for finite-length elements, and two equal moments
for zero-length elements. The Frame element considers internal deformation as
documented in the SAP2000 Analysis Reference manual. Note that P-delta effects
are not included in any constraints.

Several changes have been made to nonlinear static pushover analysis:

- P-delta effects are now included for all element types, not just for Frame
  elements, when requested in a nonlinear static pushover case. For all
  elements, P-delta effects are based on the rigid-body rotations of the
  element and do not consider the internal deformation of the element. The
  Nllink element balances the P-delta moment with a shear couple for
  finite-length elements, and two equal moments for zero-length elements. The
  Frame element does not consider internal deformation as it does for initial
  P-delta analysis. For Shell, Plane, Asolid and Solid elements, the P-delta
  effect for each step is based on the stress state at the beginning of the
  step and is not iteratively corrected at the end of the step as it is for
  Frame and Nllink elements. This means that you should specify a minimum of
  several steps for each pushover case to assure reasonable accuracy if
  P-delta effects are significant in these element types. Note that P-delta
  effects are not included in any constraints.

- Large-displacement effects are now included for all element types, not just
  for Frame elements, when requested in a nonlinear static pushover case.
  Large-displacement effects include large translations and rotations of the
  elements, but the strains within elements are still assumed to be small.
  Elements that undergo significant deformation should be divided into smaller
  elements. Note that large-displacement effects are not included in any
  constraints.

- The artificial stiffness that was being used to prevent instability when
  the stiffness of an element became zero or negative has now been reduced.
  The new procedure is more consistent, requires less iteration, is more
  likely to converge, and produces better results.

- Nonlinear static pushover analysis has been improved to reduce the
  occasional occurrence of non-convergence due to flip-flop behavior of
  yielding elements subject to P-delta effects.

- A new option is available for displacement-controlled nonlinear static
  pushover analysis. You may use the "conjugate displacement" to control the
  analysis instead of using the monitored displacement; the monitored
  displacement is still used to determine how far to push the structure and
  for plotting the pushover curve. The conjugate displacement is a weighted
  sum of all displacement degrees of freedom in the structure: each
  displacement component is multiplied by the load applied at that degree of
  freedom, and the results are summed. The conjugate displacement is usually
  the most sensitive measure of displacement in the structure under a given
  specified load. When you use the conjugate displacement to control the
  analysis, the load increments are adjusted in an attempt to reach the
  specified monitored displacement. However, the analysis will usually only
  approximately satisfy the targeted displacement, particularly if the
  monitored displacement is in a different direction than the conjugate
  displacement. The use of conjugate-displacement control is generally
  recommended.

- When using displacement control for nonlinear static pushover analyses,
  only the absolute value of the monitored displacement is used to determine
  when the analysis is complete, whether control is by the conjugate
  displacement or by the monitored displacement itself. The direction of
  loading now is determined by the specified load pattern. For example,
  positive acceleration loads will cause displacements in the negative
  direction (this is consistent with time-history analysis), regardless of
  the sign of the monitored displacement.

Several changes have been made to the Pushover Frame hinge properties:

- Separate scale factors may now be specified for the positive and negative
  axes of the stress-strain (moment-rotation or force-deformation) curves. In
  addition, the rigid-plastic nature of the hinge is shown in the plot of the
  stress-strain curves where the strains at +B and -B are subtracted out from
  the values at C, D, and E. Neither of these changes affects the results of
  previous analyses.

- For P-M-M hinges the moment-rotation curve is now interpreted as a
  relative curve, whereas in previous versions it was an absolute curve after
  scaling by the yield scale factor given for the interaction surface. Now
  when the hinge state first reaches the yield (interaction) surface, the
  moment-rotation curve is scaled by the resultant moment divided by the yield
  moment you specified at point B, i.e., it is scaled by the factor:
     (sqrt(M2*M2+M3*M3))/MB.
  Previously it was scaled by
     (sqrt(M2*M2+M3*M3))/M3max
  where M3max is the maximum value of M3 that you specified for the P-M curve
  at angle=90 deg before scaling by the yield moment. This change has no
  impact on default hinges or on hinges you defined where both the
  moment-rotation curve and interaction surface were normalized to unity at
  M3 yield. It will affect hinges where you used different scaling for the
  moment-rotation curve and the interaction surface.

- An error in the application of the moment-rotation curve in P-M-M hinges
  for points off the M3 axis has been corrected. In previous versions the
  hardening and softening tended to be under-estimated.

In version 7.21, the BS-5950 "n" factor (BS 4.3.7.6, Table 13, page 26) was
always set to one to account for destabilizing loads. This was found to be too
conservative. It has been changed back to being calculated according to Tables
13, 15, and 16. This is now the same as was being done in versions prior to
7.21 and is consistent with the SAP2000 Steel Design manual. If destabilizing
loads are present, you may overwrite the "n" factor and set it to one.

The ordering of shell joints in output from Print > Input Tables and Print >
Output Tables has been made consistent with results in the .OUT file and the
manuals. Joints 3 and 4 were previously reversed.

The on-line help has been updated to reflect the current features for
nonlinear static pushover analysis.

An error has been corrected that caused time-history reactions to neglect the
loads applied to elements that are directly connected to the reaction joint.

An error has been corrected that prevented using restraint-displacement
loading in the same model with moving load analysis.

An error has been corrected that prevented multiple analyses in the same
SAP2000 session when eigen residual-mass modes were requested.

If you prepare data using the text input data (.S2K) file, you should be
aware of the following two errors:

- In the CONSTRAINT data block the last CSYS specification was being applied
  to all constraints. This has been fixed.

- In the LOAD data block the CSYS specification is being reset to the default
  (CSYS=0) at each NAME data line. This has NOT been fixed! If you are using
  coordinate systems other than global, make sure they are specified on each
  NAME data line.

Other minor corrections and changes have been made.



Significant Changes from Version 7.11 to 7.21
---------------------------------------------

A Frameworks Plus to SAP2000 model import facility has been added. The
features currently implemented are documented in MS Word file FWP_ref.doc.
This feature is an early release and is being tested and documented before
final release.

Most input and output can now be exported in an MS Access Database format.
This feature is an early release and is being tested and documented before
final release.

The steel design under AISC-ASD and AISC-LRFD has been augmented to now also
check slender sections.

The steel design check for angle sections under AISC-ASD and AISC-LRFD has
been thoroughly revised and improved.

A concrete frame design bug that caused the shear design of columns in special
and intermediate moment resisting frames in seismic areas to be incorrect
(mostly conservative) has been corrected. This only affected the shear design
if the longitudinal reinforcement was also designed. If the column was being
checked the shear design was unaffected by this bug. Users are requested to
reverify the shear design in concrete frames if their models may be affected
by this bug.

An error has been corrected where AASHTO LRFD column shear design under
seismic load was previously not being checked for Zone 4 (severe). Users are
requested to reverify the shear design for Zone 4 if their models may be
affected by this bug.

Under BS-5950, the "n" factor (BS 4.3.7.6, Table 13, page 26) is now always
set to one to account for destabilizing loads. This is more conservative.

Other minor bug fixes have been made to the design procedures.

Longer filenames are now supported. The maximum path length is 1023
characters.

An error that sometimes caused the Printer Setup to crash on NT machines has
been fixed.

Other minor corrections and changes have been made.



Significant Changes from Version 7.10 to 7.11
---------------------------------------------

A display bug in the GUI that affected frame force diagrams when multi-valued
results were possible (Combinations, etc.) has been corrected. It only affected
the orientation of the display. All values printed and used in design were
unaffected by the bug. This bug only affected version 7.10.

A bug in AISC-LRFD steel design that affected user specified design sections
has been corrected. For most normal cases the bug made the check conservative.
This bug only affected version 7.10.

Substantial speed increases have been achieved in bridge moving-load analysis,
time-history envelope calculation, and static-pushover analysis.

Several changes have been made to the static-pushover analysis:

- An error in calculating the plastic strain for hinges in Frame elements
  with rigid-end offsets has been corrected. This error only affected elements
  with offsets at the J end of the element, and only if the offset was rigid
  (or partially so). The error could be significant in some cases. It is
  recommended that you use version 7.11 to check pushover analyses run with
  previous version of SAP2000 if your model contains rigid-end offsets.

- The static pushover analysis now checks for non-convex P-M-M interaction
  (yield) surfaces. Any P-M-M point which falls inside the surface defined
  the adjacent points is automatically moved outward (increasing the yield
  value) until a convex surface is obtained. This is an iterative process. The
  amount of error detected, and whether or not automatic correction was
  successful, is reported in the .LOG file. If significant errors are reported,
  you should correct the P-M-M surface manually.

- End releases sometimes caused excessive iteration in version 7.10. This
  has been corrected.

A minor error was corrected that may affect nonprismatic sections for Frame
elements. If the sum of the absolute lengths of all segments of a nonprismatic
section exceeds the clear length of the element, the absolute lengths are now
all scaled down so that their sum equals the clear length. This was not being
done in previous versions.



Significant Changes from Version 7.08 to 7.10
---------------------------------------------
The graphical user interface has been recompiled on a newer version of the
compiler. This should make the program run faster.

An algorithm has been added to automatically trace unbraced lengths of frame 
elements for major and minor directions of bending. The tracing accounts for 
restraints, available degrees of freedom, diaphragm constraints in XY plane, 
and support provided by other frames, shells and solids. This algorithm is 
currently used in design only for effective length K-factor calculations and 
unsupported length factors.

The limit of 2 Gigabytes on the size of the stiffness matrix has been removed
so larger models can now be analyzed. The limit is now about 16 GB.

An error was corrected in the Eurocode Steel design check where the k22, k33,
and kLTB factors were not being correctly computed. It is recommended that you
use SAP2000 version 7.10 to re-check designs from earlier version of SAP2000.

Several errors have been corrected in the AASHTO-LRFD design check. It is
recommended that you use SAP2000 version 7.10 to re-check designs from earlier
version of SAP2000.

Several enhancements have been made to the static pushover analysis from
version 7.08:

- An error in calculating the performance point introduced in version 7.08 has
  been fixed.

- Large-displacement geometric nonlinearity may now be analyzed in Frame
  elements. In particular, large rigid-body displacements and rotations are
  considered. The strains and relative bending rotations within each element
  are assumed to be small. When bending rotations are large, divide the Frame
  element into smaller elements. In addition to static pushover analysis,
  problems like snap-through buckling and cable shape-finding can be solved
  using this option. When defining a pushover case, select one of the
  following options for consideration of geometric nonlinearity effects:
  (1) None, (2) P-Delta only, or (3) P-Delta and Large displacements.

- An iterative solution strategy is now used for pushover analysis. You may now
  set the iteration tolerance and the maximum number of iterations for each
  step. The iteration tolerance is a relative force error that is allowed for
  convergence to be achieved. If convergence is not reached within the maximum
  number of iterations, the step size is halved and iteration is tried again.
  The default values are usually good enough for most analyses. Iteration is
  most important for large-displacement effects, and least important when
  geometric nonlinearities are not included in the analysis.

- You may also set the minimum number of steps, which determines the maximum
  step size allowed; the maximum number of null steps (steps of zero length
  when the solution is "turning", or steps where convergence failed), and the
  maximum number of total steps. Limiting the number of null steps prevents the
  analysis from running forever when it stalls. Limiting the number of total
  steps limits the total time of the analysis and the amount of disk space used
  to save the solution.

- For PMM hinges, you may select whether the axial force behavior is
  proportional to the specified moment-rotation curve (i.e., the entire yield
  surface expands and contracts proportionally), or is elastic-perfectly
  plastic (i.e., the yield surface does not change in the P direction). The
  first option may be more conservative and realistic. The second option may
  allow the pushover analysis to proceed further when otherwise the structure
  would fail locally under gravity load. (This feature was actually added in
  version 7.08 but was not documented in the README.TXT file)

- Minor changes have been made internally to the hinge model to increase its
  accuracy and robustness.


Significant Changes from Version 7.04 to 7.08
---------------------------------------------
SAP2000 version 7.08 provides minor enhancements to the static pushover
analysis from version 7.04. Steel design using the AISC-LRFD code now covers
slender sections and single angles.


Significant Changes from Version 7.03 to 7.04
---------------------------------------------
SAP2000 version 7.04 provides minor enhancements to the static pushover
analysis from version 7.03. Minor bugs relating to steel design have been
fixed.  


Significant Changes from Version 7.02 to 7.03
---------------------------------------------
SAP2000 version 7.03 provides minor enhancements to the static pushover
analysis from version 7.02. The original tangent approach has slightly been
modified in how it treats unloading. This version treats the unloading in a
hinge as a local redistribution. For most problems this approach should work
better in finding a solution. As in the previous approach the solution will
stall at a displacement value lower than the target if gravity load capacity
is lost in a member and other members are unable to carry this load. 
     

Significant Changes from Version 7.01 to 7.02
---------------------------------------------
SAP2000 version 7.02 provides minor enhancements and corrections to the static
pushover analysis from version 7.01, as well as a correction to the concrete
P-M-M interaction calculations for the Canadian, British and European Codes. 


Significant Changes from Version 7.00 to 7.01
---------------------------------------------
SAP2000 version 7.01 provides enhancements and corrections to the static
pushover analysis from version 7.00, as well as a correction to the moving-
load analysis. 

STATIC PUSHOVER ANALYSIS:

- Two solution strategies are now available to handle the case when a hinge
  drops load or ruptures:
  (1) The original approach that follows a tangent path to unload the 
      structure
  (2) A new secant-restart approach that goes back to the starting state
      and reloads the structure using the secant stiffness of each hinge
  The first approach is generally faster, but may fail in some cases. If that
  happens you may want to re-run the analysis using the restart approach.
  Select "Restart upon Unloading" in the Pushover Case dialog box.

- You may select whether or not you want to see negative increments that may 
  be used to unload the structure. Select "Save Positive Increments Only" in 
  the Pushover Case dialog box if you do not want to see negative increments.

- Joint displacements, frame element forces, and the hinge force-displacement 
  histories for each pushover case may be obtained from the Print Output 
  Tables option under the File menu. These are only available in spreadsheet 
  form.

- Minor corrections were made to the pushover analysis, including correction
  of the yield value for P-M-M hinges depending upon where the yield surface
  is first pierced.

- Additional information on state changes is now available in the .LOG file.

MOVING-LOAD ANALYSIS:

- An error was corrected that somtimes gave incorrect corresponding forces.
  The actual force and moment maximum and minimum values are correct.


Significant Changes from Version 6.13 to 7.00
---------------------------------------------

SAP2000 version 7.00 features a large number of new features and enhancements
as well as a few error corrections. Significant new features, enhancements and
error corrections are listed below. Minor items that do not affect the
accuracy of results are not listed.

Please see the on-line help for further information regarding all new
features and enhancements.

TOOLBAR CHANGES: The floating toolbar has been removed and replaced with a
fixed toolbar on the left-hand side of the Graphical User Interface. The main
(top) toolbar is now context-sensitive, i.e., the tool buttons available
change according to the operations available to you. See the on-line help.

DRAFTING FEATURES: The program now follows generally-accepted drafting
conventions more closely. Several new "snap" options have been added, including
snap to midpoints and ends, snap to intersections, snap to perpendicular, and
snap to lines and edges. Drawing constraints have also been added, so that
you can now draw with constrained X, Y, or Z coordinates. See the on-line help.

AERIAL VIEW: A small aerial-view window is available to enable you to see where
you are in large and complex models, and to speed up panning and zooming.
This feature may be turned on and off.

DESIGN CODES: The following codes have been added or updated to later versions:
 - Canadian CSA-A23.3-94 (Concrete)
 - AASHTO LRFD 1997 (Concrete)
 - AASHTO LRFD 1997 (Steel)
 - New Zealand NZS-3101 1995 (Concrete)
In addition, some over-conservative assumptions have been removed and the
design algorithms have been refined for all supported codes.

DXF IMPORT/EXPORT: The DXF Import and Export have been improved. You can now
import/export Solids, Joints and NLLinks and export any visible text in 2D
views only. You also now have the option to save and import data on any layer
you want. There is no longer a need to use SAP_FRAMES and SAP_SHELLS. You must
use a template DXF file when exporting. The template file can include layer
names, color settings etc. There is a template file provided. The on-line help
has the full documentation of the new DXF import and Export.

WINDOWS METAFILE EXPORT: You may now export graphics in the Windows Enhanced
Metafile format for inclusion in documents using other programs that support
this format.

COLOR PRINTER SETTINGS: There are new color-printer settings available in the
Options... Color Setting menu item.

ORTHOTROPIC MATERIAL: Orthotropic materials may now be defined and modified
in the Graphical User Interface (they have always been available through the
text file). See the on-line help.

ORTHOTROPIC SHELL: The Shell element now supports orthotropic material
properties. Orthotropic properties may be rotated through a "material angle"
when defining Shell Section properties. See the on-line help.

THICK PLATE: A thick-plate formulation that includes transverse-shear
deformation is now available for the Shell element. The previous thin-plate
formulation that neglects transverse-shear deformation is still available. You
can select whether or not to use the thick-plate formulation in the Sehll
Section properties definition. Only plate behavior is affected, not membrane
behavior. See the on-line help. THE THICK-PLATE OPTION IS NOW THE DEFAULT.

The thick-plate formulation tends to be somewhat stiffer yet more accurate
than the thin-plate formulation for most thin-plate bending problems. Of
course, the presence of shearing deformations will make the thick-plate
formulation more flexible.

RESIDUAL-MASS MODES: Residual-mass (missing-mass) modes may now be calculated
as an option for eigen-analysis. When present, residual-mass modes are
automatically included in Response-spectrum and Time-history analyses. For Ritz
analysis, residual mass has always been included for all starting load vectors.
This has not been changed. See the on-line help.

EFFECTIVE-PERIOD PRINTOUT: For each load case, the effective period is
printed in the output (.OUT) file under the heading MODAL LOAD PARTICIPATION
RATIOS. The effective period is calculated as the Rayleigh period of the
structure using the static deflection for the load case as the mode shape. This
information can be used for determining an equivalent single-degree-of-freedom
system for blast-load analysis and other special purposes. It is not of
general utility.

STATIC PUSHOVER ANALYSIS: Static Pushover analysis has been added to the
Nonlinear version of SAP2000. Nonlinear hinges may be defined anywhere in the
Frame elements. Properties may be user-defined or may be calculated
automatically by the program from Section Properties according to standard
pushover guidelines. Analyses may be force or displacement controlled. Pushover
results may be displayed graphically and in table form. See the on-line help.

KEY-DEVICE DEPARTMENT NAME: For large organizations working over a wide-area
network, you may define "department names" to restrict the use of the hardware
key device to certain people. Please contact CSI for more information if you
need to use this feature.

TIME-HISTORY ELEMENT LOADS: Time-history results now include the effect of
all internal element loads, including gravity, concentrated, distributed and
prestress loads on Frame elements, as well as temperature loads on all element
types.

BRIDGE P-DELTA FORCES: An error was corrected that sometimes gave incorrect
Frame element forces for Moving-load cases when P-Delta forces are present.


Significant Changes from Version 6.11 to 6.13
---------------------------------------------

GENERAL...

(a) An error was corrected that produced incorrect moments in Frame
    elements subjected to trapezoidal loads that varied within an output
    segment. The error was zero if changes in the load corresponded to
    output locations.

(b) An error was corrected that calculated the mass incorrectly for
    non-prismatic Frame elements. This error was generally small unless
    the length of the element was less than unity in the length units
    used for the analysis.

(c) An error was corrected that sometimes produced incorrect results for
    models that have both Welds and Constraints.

(d) An error was corrected in which gravity loads imported from SAP90
    text files were not rotated to account for changing global up to +Z.
    This does not affect SAP90 files for which global +Z is already
    upward.

DESIGN...

(e) An error was corrected where the program was not accounting for the
    optional user-specified flexural steel provided when calculating
    flexural capacity for shear design of beams where ductile design is
    requested. Version 6.11 shear design was based on required values
    rather than the user-provided flexural values.

(f) An error was corrected for concrete design by the Canadian code
    where the minimum beam reinforcement on the tension side was
    unconditionally related to (4/3)As. This has been corrected for
    seismic design and is now consistent with the SAP2000 Concrete
    Design Manual.

(g) An error was corrected where 1% minimum steel required by the ACI
    code for concrete column design was being used for all codes, which
    is over-conservative for some codes. This has been corrected to be
    code-dependent.

(h) For ductile design, the required minimum in-span moment capacities
    for the ACI and Canadian codes has been slightly refined. This
    should have minimal impact on previous designs.

BRIDGE ANALYSIS...

(i) An error was corrected that could produce incorrect moving-load
    results for lanes that contain Frame elements with their local 1
    axes oriented in the opposite direction of the direction of the
    lane.

    For example, suppose element 1 connects to joints 1 and 2, and
    element 2 connects to joints 2 and 3. Results in version 6.11 for a
    Lane defined from element 1 to element 2 would have been correct,
    but for a Lane defined from element 2 to element 1 may have been
    incorrect. Results in version 6.13 should be correct for both cases.

(j) An error was corrected that sometimes produced incorrect
    corresponding quantities in a moving-load analysis. For example, the
    maximum value of M3 would be correct, but the value of M3
    corresponding to maximum M2 may have been incorrect.

(k) Moving load analysis has been speeded up significantly for some
    models involving multiple lanes and vehicles.

(l) An error was corrected that sometimes produced an incorrect
    display of an influence line.

TIME-HISTORY ANALYSIS...

(m) An error was corrected that caused incorrect importing from SAP90
    and SAP2000 text files of differing modal damping values for a
    time-history analysis. The error did not affect constant modal
    damping for all modes in a time-history analysis, nor differing
    modal damping values specified within the graphical user interface.

NONLINEAR TIME-HISTORY ANALYSIS...

(o) An error was corrected that sometimes displayed Nllink time-history traces
    for the wrong Nllink element.

(p) An error was corrected that sometimes produced incorrect History envelopes
    for Nllink elements.

(q) An error was corrected that produced incorrect reactions for nonlinear
    time-history analyses at joints connected to Nllink elements. Envelopes and
    time-history traces of the reactions were both calculated using the linear
    effective stiffness rather than the nonlinear element properties. Only
    joints connected to Nllinks were affected, and only during nonlinear
    time-history analyses.

(r) An error was corrected that produced incorrect group force sums for
    nonlinear time-history analyses for groups containing Nllink elements.

MISCELLANEOUS...

(s) Minor improvements and corrections were made to the Graphical User
    Interface to fix problems with the functioning of the dialog boxes
    and the display of results, and to improve the stability of the
    program. These changes should not have any significant effect on
    previous analysis or design results.




=====================================================
6. Upgrading to SAP2000 from SAP90 Version 5.4 or 5.5
=====================================================

Most modeling and analysis features available in SAP90 are also present in
SAP2000, and many new features have been added. Only the SAP90 heat-transfer
analysis features are not currently available in SAP2000.

SAP90 input data files for versions 5.4 and 5.5 may be imported directly into
SAP2000. A translated SAP2000 input data text file (with extension .S2K) will
be created, and the model will be stored in the native SAP2000 data base file
(with extension .SDB). The model can then be modified, analyzed, designed, and
displayed.

When you import a SAP90 input data file, SAP2000 will ask you to specify what
direction was assumed to be upward in the SAP90 model. All coordinate-dependent
quantities in the SAP90 model will be converted to conform with the SAP2000
convention that the +Z direction is upward.

Please note that SAP90 Bridge Moving Load cases may be translated into
multiple SAP2000 Moving Load cases. In many cases, you may be able to combine
them back into a single SAP2000 Moving Load case, but you will have to do this
manually. For each SAP90 Moving Load case, a single SAP2000 envelope- type
Combo is created during translation. Please also note that the "pm" portion of
the SAP90 lane load was used for moments M2 and M3 for all elements, while the
floating concentrated load "pm" for the SAP2000 vehicle load only applies to
the vertical moment (usually M3), and only for lane elements!

WARNING! Some imported data may be interpreted differently by SAP2000 than by
SAP90. For example, the interaction between end offsets and end releases is
different, as is the interaction between prestress load and P-Delta analysis.
Be sure to check your imported model carefully! Compare the results of
analyses using both SAP90 and SAP2000 before making further use of the
imported model!

See Topic "Upgrading from SAP90" in Chapter II of the "SAP2000 Getting
Started" manual for more information.




=======================
7. Analysis Executables
=======================

With SAP2000 you can create your model, perform the analysis, display the
results, and check the design all within a single graphical user interface.
However, some users have requested the ability to perform multiple analyses in
a batch mode outside of the user interface. We have provided executable files
for this purpose in a subfolder called UTILITY under the SAP2000 folder.

These files are 32-bit Windows executables. They can be run from a Command
Prompt (DOS-like) window under Windows 95/98/NT/2000. These are not DOS
executables and cannot be run on DOS or Windows 3.x machines. In order to run
properly, they must be copied or moved from the UTILITY subfolder up to the
SAP2000 folder where SAP2000.EXE is located.

The executables are:

(a) SAPRE.EXE. This reads and checks a SAP2000 input text file and prepares
    the data for subsequent analysis using SAPGO.EXE. To run, give the
    command:

       SAPRE  sap2000file  /M:nnnnn

    where "sap2000file" (required) is the SAP2000 input text file to be
    read; and "nnnnn" (optional) is the amount of memory to be allocated
    in 1000-byte units. The "sap2000file" must include the extension .S2K.
    The default memory if the /M: parameter is omitted is the same as
    specifying /M:2000.

(b) SAPGO.EXE. This perform the analysis following the execution of SAPRE.EXE.
    To run, give the command:

       SAPGO  sap2000file  /M:nnnnn

    where "sap2000file" (required) is the SAP2000 input text file to be
    analyzed; and "nnnnn" (optional) is the amount of memory to be allocated
    in 1000-byte units. The default memory if the /M: parameter is omitted is
    the same as specifying 1/16 of physical memory.

To perform an analysis using these executables:

(a) Prepare an input text file (.S2K) using a text editor or by exporting a
    file from the SAP2000 graphical user interface.

(b) Run SAPRE.EXE followed by SAPGO.EXE. These commands may be place in a
    batch (.BAT) file.

(c) To display the results, import the SAP2000 .JOB file under the File menu
    of the SAP2000 graphical user interface. It is important to import the
    .JOB file, not the .S2K file, or your analysis results will be lost.
