#############################################################################
# Generated by PAGE version 4.18
#  in conjunction with Tcl version 8.6
#  Jan 16, 2019 11:11:27 PM CET  platform: Linux
set vTcl(timestamp) ""


if {!$vTcl(borrow)} {

vTcl:font:add_GUI_font TkDefaultFont \
"-family Helvetica -size -12 -weight normal -slant roman -underline 0 -overstrike 0"
vTcl:font:add_GUI_font TkFixedFont \
"-family courier -size -12 -weight normal -slant roman -underline 0 -overstrike 0"
vTcl:font:add_GUI_font TkMenuFont \
"-family Helvetica -size -12 -weight normal -slant roman -underline 0 -overstrike 0"
vTcl:font:add_GUI_font TkTextFont \
"-family Helvetica -size -12 -weight normal -slant roman -underline 0 -overstrike 0"
set vTcl(actual_gui_font_dft_name) TkDefaultFont
set vTcl(actual_gui_font_menu_name) TkMenuFont
set vTcl(actual_gui_bg) #d9d9d9
set vTcl(actual_gui_fg) #000000
set vTcl(actual_gui_menu_bg) #d9d9d9
set vTcl(actual_gui_menu_fg) #000000
set vTcl(complement_color) #d9d9d9
set vTcl(analog_color_p) #d9d9d9
set vTcl(analog_color_m) #d9d9d9
set vTcl(active_fg) #000000
set vTcl(actual_gui_menu_active_bg)  #d8d8d8
set vTcl(active_menu_fg) #000000
}

#############################################################################
# vTcl Code to Load User Fonts

vTcl:font:add_font \
    "-family gothic -size 9 -weight normal -slant roman -underline 0 -overstrike 0" \
    user \
    vTcl:font10
vTcl:font:add_font \
    "-family fixed -size 10 -weight bold -slant roman -underline 0 -overstrike 0" \
    user \
    vTcl:font9
#################################
#LIBRARY PROCEDURES
#


if {[info exists vTcl(sourcing)]} {

proc vTcl:project:info {} {
    set base .top39
    global vTcl
    set base $vTcl(btop)
    if {$base == ""} {
        set base .top39
    }
    namespace eval ::widgets::$base {
        set dflt,origin 0
        set runvisible 1
    }
    namespace eval ::widgets_bindings {
        set tagslist _TopLevel
    }
    namespace eval ::vTcl::modules::main {
        set procs {
        }
        set compounds {
        }
        set projectType single
    }
}
}

#################################
# GENERATED GUI PROCEDURES
#

proc vTclWindow.top39 {base} {
    if {$base == ""} {
        set base .top39
    }
    if {[winfo exists $base]} {
        wm deiconify $base; return
    }
    set top $base
    ###################
    # CREATING WIDGETS
    ###################
    vTcl::widgets::core::toplevel::createCmd $top -class Toplevel \
        -relief sunken -background {#d9d9d9} -highlightcolor black 
    wm focusmodel $top passive
    wm geometry $top 601x728+472+131
    update
    # set in toplevel.wgt.
    global vTcl
    global img_list
    set vTcl(save,dflt,origin) 0
    wm maxsize $top 1905 1050
    wm minsize $top 1 1
    wm overrideredirect $top 0
    wm resizable $top 1 1
    wm deiconify $top
    wm title $top "Pyxtal Main Controls"
    vTcl:DefineAlias "$top" "pyxtal" vTcl:Toplevel:WidgetProc "" 1
    labelframe $top.lab41 \
        -foreground black -text {Input Files} -background {#d9d9d9} \
        -height 335 -highlightcolor black -width 290 
    vTcl:DefineAlias "$top.lab41" "InputFileFrame" vTcl:WidgetProc "pyxtal" 1
    set site_3_0 $top.lab41
    vTcl::widgets::ttk::scrolledlistbox::CreateCmd $site_3_0.scr42 \
        -background {#d9d9d9} -height 75 -highlightcolor black -width 125 
    vTcl:DefineAlias "$site_3_0.scr42" "fileListbox" vTcl:WidgetProc "pyxtal" 1

    $site_3_0.scr42.01 configure -background white \
        -font TkFixedFont \
        -foreground black \
        -height 3 \
        -highlightcolor #d9d9d9 \
        -selectbackground #c4c4c4 \
        -selectforeground black \
        -width 10
    button $site_3_0.but93 \
        -activebackground {#d9d9d9} -activeforeground black \
        -background {#d9d9d9} -command addButtonCommand -foreground {#000000} \
        -highlightcolor black -text {Add Files} 
    vTcl:DefineAlias "$site_3_0.but93" "addButton" vTcl:WidgetProc "pyxtal" 1
    button $site_3_0.but94 \
        -activebackground {#d9d9d9} -activeforeground black \
        -background {#d9d9d9} -command clearButtonCommand \
        -foreground {#000000} -highlightcolor black -text {Clear Files} 
    vTcl:DefineAlias "$site_3_0.but94" "clearButton" vTcl:WidgetProc "pyxtal" 1
    text $site_3_0.tex40 \
        -background white -font TkTextFont -foreground black -height 30 \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black -width 266 \
        -wrap word 
    .top39.lab41.tex40 configure -font TkTextFont
    .top39.lab41.tex40 insert end text
    vTcl:DefineAlias "$site_3_0.tex40" "pathBox" vTcl:WidgetProc "pyxtal" 1
    label $site_3_0.lab41 \
        -activebackground {#f9f9f9} -activeforeground black \
        -background {#d9d9d9} -font $::vTcl(fonts,vTcl:font9,object) \
        -foreground {#000000} -highlightcolor black -text Files: 
    vTcl:DefineAlias "$site_3_0.lab41" "filesLabel" vTcl:WidgetProc "pyxtal" 1
    label $site_3_0.lab42 \
        -activebackground {#f9f9f9} -activeforeground black \
        -background {#d9d9d9} -font $::vTcl(fonts,vTcl:font9,object) \
        -foreground {#000000} -highlightcolor black -text Path: 
    vTcl:DefineAlias "$site_3_0.lab42" "pathLabel" vTcl:WidgetProc "pyxtal" 1
    place $site_3_0.scr42 \
        -in $site_3_0 -x 10 -y 160 -width 276 -relwidth 0 -height 170 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.but93 \
        -in $site_3_0 -x 10 -y 30 -width 98 -height 33 -anchor nw \
        -bordermode ignore 
    place $site_3_0.but94 \
        -in $site_3_0 -x 170 -y 30 -width 108 -height 33 -anchor nw \
        -bordermode ignore 
    place $site_3_0.tex40 \
        -in $site_3_0 -x 10 -y 100 -width 266 -relwidth 0 -height 30 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.lab41 \
        -in $site_3_0 -x 10 -y 140 -anchor nw -bordermode ignore 
    place $site_3_0.lab42 \
        -in $site_3_0 -x 10 -y 80 -width 46 -height 17 -anchor nw \
        -bordermode ignore 
    button $top.but43 \
        -activebackground {#d9d9d9} -activeforeground black \
        -background {#d9d9d9} -command saveButtonCommand \
        -foreground {#000000} -highlightcolor black -text {Save Defaults} 
    vTcl:DefineAlias "$top.but43" "saveDefButton" vTcl:WidgetProc "pyxtal" 1
    labelframe $top.lab48 \
        -foreground black -text {Input File Options} -background {#d9d9d9} \
        -height 285 -highlightcolor black -width 290 
    vTcl:DefineAlias "$top.lab48" "inputOptionFrame" vTcl:WidgetProc "pyxtal" 1
    set site_3_0 $top.lab48
    radiobutton $site_3_0.rad41 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -command inFileTypeChange -foreground {#000000} \
        -highlightcolor black -justify left -text {Image file} -value image \
        -variable {} 
    vTcl:DefineAlias "$site_3_0.rad41" "imageRadio" vTcl:WidgetProc "pyxtal" 1
    radiobutton $site_3_0.rad42 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -command inFileTypeChange -foreground {#000000} \
        -highlightcolor black -justify left -text {gsd assemblies} \
        -value assemblies -variable {} 
    vTcl:DefineAlias "$site_3_0.rad42" "gsdAssemRadio" vTcl:WidgetProc "pyxtal" 1
    radiobutton $site_3_0.rad43 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -command inFileTypeChange -foreground {#000000} \
        -highlightcolor black -justify left -text {gsd particles} \
        -value particles -variable {} 
    vTcl:DefineAlias "$site_3_0.rad43" "gsdPartRadio" vTcl:WidgetProc "pyxtal" 1
    entry $site_3_0.ent44 \
        -background white -font TkFixedFont -foreground {#000000} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable sphereSizeStr 
    vTcl:DefineAlias "$site_3_0.ent44" "sphereEntry" vTcl:WidgetProc "pyxtal" 1
    bind $site_3_0.ent44 <FocusOut> {
        lambda e:pyxtal_support.validateInteger(e,
                                    self.sphereSizeStr,
                                    self.sphereSize)
    }
    label $site_3_0.lab45 \
        -activebackground {#f9f9f9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -text {Sphere Size:} 
    vTcl:DefineAlias "$site_3_0.lab45" "sphereSizeLabel" vTcl:WidgetProc "pyxtal" 1
    label $site_3_0.lab53 \
        -activebackground {#f9f9f9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -text {of type:} 
    vTcl:DefineAlias "$site_3_0.lab53" "partTypeLabel" vTcl:WidgetProc "pyxtal" 1
    entry $site_3_0.ent54 \
        -background white -font TkFixedFont -foreground {#000000} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable partTypeStr 
    vTcl:DefineAlias "$site_3_0.ent54" "partTypeEntry" vTcl:WidgetProc "pyxtal" 1
    checkbutton $site_3_0.che67 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -offrelief sunken -text {Periodic Boundaries} \
        -variable periodBound 
    vTcl:DefineAlias "$site_3_0.che67" "periodicCheck" vTcl:WidgetProc "pyxtal" 1
    labelframe $site_3_0.lab82 \
        -foreground black -text frames -background {#d9d9d9} -height 75 \
        -highlightcolor black -width 270 
    vTcl:DefineAlias "$site_3_0.lab82" "framesFrame" vTcl:WidgetProc "pyxtal" 1
    set site_4_0 $site_3_0.lab82
    entry $site_4_0.ent83 \
        -background white -font TkFixedFont -foreground {#000000} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable fromFrameStr 
    vTcl:DefineAlias "$site_4_0.ent83" "fromEntry" vTcl:WidgetProc "pyxtal" 1
    bind $site_4_0.ent83 <FocusOut> {
        lambda e:pyxtal_support.validateInteger(e,
                                    self.fromFrameStr,
                                    self.fromFrame)
    }
    label $site_4_0.lab84 \
        -activebackground {#f9f9f9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -text From 
    vTcl:DefineAlias "$site_4_0.lab84" "fromLabel" vTcl:WidgetProc "pyxtal" 1
    label $site_4_0.lab85 \
        -activebackground {#f9f9f9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -text to 
    vTcl:DefineAlias "$site_4_0.lab85" "toLabel" vTcl:WidgetProc "pyxtal" 1
    entry $site_4_0.ent86 \
        -background white -font TkFixedFont -foreground {#000000} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable toFrameStr 
    vTcl:DefineAlias "$site_4_0.ent86" "toEntry" vTcl:WidgetProc "pyxtal" 1
    bind $site_4_0.ent86 <FocusOut> {
        lambda e:pyxtal_support.validateInteger(e,
                                    self.toFrameStr,
                                    self.toFrame)
    }
    label $site_4_0.lab87 \
        -activebackground {#f9f9f9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -text by 
    vTcl:DefineAlias "$site_4_0.lab87" "byLabel" vTcl:WidgetProc "pyxtal" 1
    entry $site_4_0.ent88 \
        -background white -font TkFixedFont -foreground {#000000} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable byFrameStr 
    vTcl:DefineAlias "$site_4_0.ent88" "byEntry" vTcl:WidgetProc "pyxtal" 1
    bind $site_4_0.ent88 <FocusOut> {
        lambda e:pyxtal_support.validateInteger(e,
                                    self.byFrameStr,
                                    self.byFrame)
    }
    place $site_4_0.ent83 \
        -in $site_4_0 -x 50 -y 35 -width 36 -height 27 -anchor nw \
        -bordermode ignore 
    place $site_4_0.lab84 \
        -in $site_4_0 -x 10 -y 40 -width 38 -relwidth 0 -height 15 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_4_0.lab85 \
        -in $site_4_0 -x 100 -y 40 -width 38 -height 15 -anchor nw \
        -bordermode ignore 
    place $site_4_0.ent86 \
        -in $site_4_0 -x 130 -y 35 -width 36 -height 27 -anchor nw \
        -bordermode ignore 
    place $site_4_0.lab87 \
        -in $site_4_0 -x 190 -y 40 -width 38 -height 15 -anchor nw \
        -bordermode ignore 
    place $site_4_0.ent88 \
        -in $site_4_0 -x 220 -y 35 -width 36 -height 27 -anchor nw \
        -bordermode ignore 
    checkbutton $site_3_0.che92 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text {Dark Spheres} -variable darkSpheres 
    vTcl:DefineAlias "$site_3_0.che92" "darkSpheresCheck" vTcl:WidgetProc "pyxtal" 1
    place $site_3_0.rad41 \
        -in $site_3_0 -x 10 -y 30 -width 144 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.rad42 \
        -in $site_3_0 -x 10 -y 90 -width 159 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.rad43 \
        -in $site_3_0 -x 10 -y 60 -width 132 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.ent44 \
        -in $site_3_0 -x 120 -y 215 -width 66 -relwidth 0 -height 27 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.lab45 \
        -in $site_3_0 -x 15 -y 220 -width 98 -relwidth 0 -height 15 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.lab53 \
        -in $site_3_0 -x 150 -y 90 -width 68 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.ent54 \
        -in $site_3_0 -x 210 -y 86 -width 66 -relwidth 0 -height 27 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.che67 \
        -in $site_3_0 -x 10 -y 250 -width 185 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.lab82 \
        -in $site_3_0 -x 10 -y 120 -width 270 -relwidth 0 -height 75 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.che92 \
        -in $site_3_0 -x 150 -y 30 -width 135 -height 17 -anchor nw \
        -bordermode ignore 
    button $top.but40 \
        -activebackground {#d9d9d9} -activeforeground black \
        -background {#d9d9d9} -command loadButtonCommand \
        -foreground {#000000} -highlightcolor black -text {Load Defaults} 
    vTcl:DefineAlias "$top.but40" "loadDefButton" vTcl:WidgetProc "pyxtal" 1
    labelframe $top.lab46 \
        -foreground black -text {Output File Options} -background {#d9d9d9} \
        -height 185 -highlightcolor black -width 270 
    vTcl:DefineAlias "$top.lab46" "outputFrame" vTcl:WidgetProc "pyxtal" 1
    set site_3_0 $top.lab46
    entry $site_3_0.ent44 \
        -background white -font TkFixedFont -foreground {#000000} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable imageSizeStr 
    vTcl:DefineAlias "$site_3_0.ent44" "imageSizeEntry" vTcl:WidgetProc "pyxtal" 1
    bind $site_3_0.ent44 <FocusOut> {
        lambda e:pyxtal_support.validateInteger(e,
                                    self.imageSizeStr,
                                    self.imageSize)
    }
    label $site_3_0.lab47 \
        -activebackground {#f9f9f9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -text {Output Image size:} 
    vTcl:DefineAlias "$site_3_0.lab47" "imageSizeLabel" vTcl:WidgetProc "pyxtal" 1
    checkbutton $site_3_0.che68 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text {Image and Circles} -variable outCircles 
    vTcl:DefineAlias "$site_3_0.che68" "outCirclesCheck" vTcl:WidgetProc "pyxtal" 1
    checkbutton $site_3_0.che69 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text Triangulation -variable outTriang 
    vTcl:DefineAlias "$site_3_0.che69" "outTriangCheck" vTcl:WidgetProc "pyxtal" 1
    checkbutton $site_3_0.che70 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text {Image + angle + defects} -variable outAll 
    vTcl:DefineAlias "$site_3_0.che70" "outAllCheck" vTcl:WidgetProc "pyxtal" 1
    checkbutton $site_3_0.che91 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text mpeg -variable outMpeg 
    vTcl:DefineAlias "$site_3_0.che91" "outMpegCheck" vTcl:WidgetProc "pyxtal" 1
    place $site_3_0.ent44 \
        -in $site_3_0 -x 160 -y 110 -width 76 -relwidth 0 -height 27 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.lab47 \
        -in $site_3_0 -x 15 -y 120 -width 128 -relwidth 0 -height 15 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.che68 \
        -in $site_3_0 -x 10 -y 30 -width 175 -height 17 -anchor nw \
        -bordermode ignore 
    place $site_3_0.che69 \
        -in $site_3_0 -x 10 -y 60 -width 175 -height 17 -anchor nw \
        -bordermode ignore 
    place $site_3_0.che70 \
        -in $site_3_0 -x 10 -y 90 -width 205 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.che91 \
        -in $site_3_0 -x 10 -y 150 -width 105 -height 17 -anchor nw \
        -bordermode ignore 
    labelframe $top.lab57 \
        -foreground black -text Analysis -background {#d9d9d9} -height 245 \
        -highlightcolor black -width 270 
    vTcl:DefineAlias "$top.lab57" "analysisFrame" vTcl:WidgetProc "pyxtal" 1
    set site_3_0 $top.lab57
    checkbutton $site_3_0.che64 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text {Sphere Trajectories} -variable doTraject 
    vTcl:DefineAlias "$site_3_0.che64" "trajCheck" vTcl:WidgetProc "pyxtal" 1
    checkbutton $site_3_0.che71 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text {Orientational Correlation func.} \
        -variable doOrientCorr 
    vTcl:DefineAlias "$site_3_0.che71" "orientCorrCheck" vTcl:WidgetProc "pyxtal" 1
    checkbutton $site_3_0.che43 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text {Area Orientation Histogram} \
        -variable doOrientHist 
    vTcl:DefineAlias "$site_3_0.che43" "orientHistCheck" vTcl:WidgetProc "pyxtal" 1
    checkbutton $site_3_0.che44 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text {Sphere Statistics} -variable doSphereStats 
    vTcl:DefineAlias "$site_3_0.che44" "sphereStatsCheck" vTcl:WidgetProc "pyxtal" 1
    checkbutton $site_3_0.che45 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text {Chrystal Defect Statistics} \
        -variable doDefectStats 
    vTcl:DefineAlias "$site_3_0.che45" "defectStatsCheck" vTcl:WidgetProc "pyxtal" 1
    checkbutton $site_3_0.che46 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text {Depth Profile (z axis)} -variable doZProfile 
    vTcl:DefineAlias "$site_3_0.che46" "zProfileCheck" vTcl:WidgetProc "pyxtal" 1
    checkbutton $site_3_0.che47 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text {Meaning of Life} -variable doMeaningLife 
    vTcl:DefineAlias "$site_3_0.che47" "meaningLifeCheck" vTcl:WidgetProc "pyxtal" 1
    place $site_3_0.che64 \
        -in $site_3_0 -x 10 -y 90 -width 205 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.che71 \
        -in $site_3_0 -x 10 -y 180 -width 245 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.che43 \
        -in $site_3_0 -x 10 -y 150 -width 225 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.che44 \
        -in $site_3_0 -x 10 -y 60 -width 195 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.che45 \
        -in $site_3_0 -x 10 -y 120 -width 215 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.che46 \
        -in $site_3_0 -x 10 -y 30 -width 195 -height 17 -anchor nw \
        -bordermode ignore 
    place $site_3_0.che47 \
        -in $site_3_0 -x 10 -y 210 -width 245 -height 17 -anchor nw \
        -bordermode ignore 
    labelframe $top.lab97 \
        -foreground black -text {Window Control} -background {#d9d9d9} \
        -height 145 -highlightcolor black -width 270 
    vTcl:DefineAlias "$top.lab97" "windowFrame" vTcl:WidgetProc "pyxtal" 1
    set site_3_0 $top.lab97
    checkbutton $site_3_0.che98 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text {Retain Windows} -variable retainWin 
    vTcl:DefineAlias "$site_3_0.che98" "retainCheck" vTcl:WidgetProc "pyxtal" 1
    button $site_3_0.but99 \
        -activebackground {#d9d9d9} -activeforeground black \
        -background {#d9d9d9} -command GoButtonCommand \
        -font $::vTcl(fonts,vTcl:font10,object) -foreground {#000000} \
        -highlightcolor black -text Go 
    vTcl:DefineAlias "$site_3_0.but99" "goButton" vTcl:WidgetProc "pyxtal" 1
    button $site_3_0.but100 \
        -activebackground {#d9d9d9} -activeforeground black \
        -background {#d9d9d9} -command killAllButtonCommand \
        -foreground {#000000} -highlightcolor black \
        -text {Kill All 
Image Windows} 
    vTcl:DefineAlias "$site_3_0.but100" "killAllButton" vTcl:WidgetProc "pyxtal" 1
    checkbutton $site_3_0.che56 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text {Lock zoom} -variable lockZoom 
    vTcl:DefineAlias "$site_3_0.che56" "lockZoomCheck" vTcl:WidgetProc "pyxtal" 1
    checkbutton $site_3_0.che57 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text {Lock views} -variable lockViews 
    vTcl:DefineAlias "$site_3_0.che57" "lockViewsCheck" vTcl:WidgetProc "pyxtal" 1
    checkbutton $site_3_0.che43 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -foreground {#000000} -highlightcolor black \
        -justify left -text {Quiet Batchmode} -variable batchmode -width 129 
    vTcl:DefineAlias "$site_3_0.che43" "batchmodeCheck" vTcl:WidgetProc "pyxtal" 1
    place $site_3_0.che98 \
        -in $site_3_0 -x 10 -y 60 -width 119 -height 17 -anchor nw \
        -bordermode ignore 
    place $site_3_0.but99 \
        -in $site_3_0 -x 10 -y 90 -width 108 -relwidth 0 -height 43 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.but100 \
        -in $site_3_0 -x 150 -y 90 -width 108 -relwidth 0 -height 43 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.che56 \
        -in $site_3_0 -x 150 -y 60 -width 109 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.che57 \
        -in $site_3_0 -x 150 -y 30 -width 109 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.che43 \
        -in $site_3_0 -x 10 -y 30 -width 129 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    labelframe $top.lab50 \
        -foreground black -text Progress -background {#d9d9d9} -height 75 \
        -highlightcolor black -width 270 
    vTcl:DefineAlias "$top.lab50" "progressFrame" vTcl:WidgetProc "pyxtal" 1
    set site_3_0 $top.lab50
    ttk::progressbar $site_3_0.tPr51
    vTcl:DefineAlias "$site_3_0.tPr51" "fileProgressbar" vTcl:WidgetProc "pyxtal" 1
    message $site_3_0.mes52 \
        -anchor w -background {#d9d9d9} -foreground {#000000} \
        -highlightcolor black -text {Files: } -width 131 
    vTcl:DefineAlias "$site_3_0.mes52" "fileMessage" vTcl:WidgetProc "pyxtal" 1
    message $site_3_0.mes53 \
        -anchor w -background {#d9d9d9} -foreground {#000000} \
        -highlightcolor black -text {Frames: } -width 151 
    vTcl:DefineAlias "$site_3_0.mes53" "frameMessage" vTcl:WidgetProc "pyxtal" 1
    ttk::progressbar $site_3_0.tPr54
    vTcl:DefineAlias "$site_3_0.tPr54" "frameProgressbar" vTcl:WidgetProc "pyxtal" 1
    place $site_3_0.tPr51 \
        -in $site_3_0 -x 150 -y 20 -width 110 -relwidth 0 -height 19 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.mes52 \
        -in $site_3_0 -x 10 -y 20 -width 131 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $site_3_0.mes53 \
        -in $site_3_0 -x 10 -y 50 -width 151 -height 17 -anchor nw \
        -bordermode ignore 
    place $site_3_0.tPr54 \
        -in $site_3_0 -x 150 -y 50 -width 110 -relwidth 0 -height 19 \
        -relheight 0 -anchor nw -bordermode ignore 
    ###################
    # SETTING GEOMETRY
    ###################
    place $top.lab41 \
        -in $top -x 10 -y 80 -width 290 -relwidth 0 -height 335 -relheight 0 \
        -anchor nw -bordermode ignore 
    place $top.but43 \
        -in $top -x 170 -y 30 -width 108 -relwidth 0 -height 33 -relheight 0 \
        -anchor nw -bordermode ignore 
    place $top.lab48 \
        -in $top -x 10 -y 430 -width 290 -relwidth 0 -height 285 -relheight 0 \
        -anchor nw -bordermode ignore 
    place $top.but40 \
        -in $top -x 30 -y 30 -width 108 -relwidth 0 -height 33 -relheight 0 \
        -anchor nw -bordermode ignore 
    place $top.lab46 \
        -in $top -x 320 -y 20 -width 270 -relwidth 0 -height 185 -relheight 0 \
        -anchor nw -bordermode ignore 
    place $top.lab57 \
        -in $top -x 320 -y 220 -width 270 -relwidth 0 -height 245 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $top.lab97 \
        -in $top -x 320 -y 480 -width 270 -relwidth 0 -height 145 \
        -relheight 0 -anchor nw -bordermode ignore 
    place $top.lab50 \
        -in $top -x 320 -y 640 -width 270 -relwidth 0 -height 75 -relheight 0 \
        -anchor nw -bordermode ignore 

    vTcl:FireEvent $base <<Ready>>
}

#############################################################################
## Binding tag:  _TopLevel

bind "_TopLevel" <<Create>> {
    if {![info exists _topcount]} {set _topcount 0}; incr _topcount
}
bind "_TopLevel" <<DeleteWindow>> {
    if {[set ::%W::_modal]} {
                vTcl:Toplevel:WidgetProc %W endmodal
            } else {
                destroy %W; if {$_topcount == 0} {exit}
            }
}
bind "_TopLevel" <Destroy> {
    if {[winfo toplevel %W] == "%W"} {incr _topcount -1}
}

set btop ""
if {$vTcl(borrow)} {
    set btop .bor[expr int([expr rand() * 100])]
    while {[lsearch $btop $vTcl(tops)] != -1} {
        set btop .bor[expr int([expr rand() * 100])]
    }
}
set vTcl(btop) $btop
Window show .
Window show .top39 $btop
if {$vTcl(borrow)} {
    $btop configure -background plum
}

