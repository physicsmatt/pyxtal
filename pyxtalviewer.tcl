#############################################################################
# Generated by PAGE version 4.17
# in conjunction with Tcl version 8.6
# Oct 19, 2018 11:13:48 AM CEST  platform: Linux
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
        -background {#d9d9d9} -highlightcolor black 
    wm focusmodel $top passive
    wm geometry $top 822x867+549+122
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
    wm title $top "Pyxtal Viewer"
    vTcl:DefineAlias "$top" "Toplevel1" vTcl:Toplevel:WidgetProc "" 1
    labelframe $top.lab40 \
        -foreground black -text Image -background {#d9d9d9} -height 45 \
        -highlightcolor black -width 280 
    vTcl:DefineAlias "$top.lab40" "imageframe" vTcl:WidgetProc "Toplevel1" 1
    set site_3_0 $top.lab40
    radiobutton $site_3_0.rad45 \
        -activebackground {#d9d9d9} -activeforeground black \
        -background {#d9d9d9} -command showImageChange -foreground {#000000} \
        -highlightcolor black -justify left -text Raw -value raw \
        -variable whichImage 
    vTcl:DefineAlias "$site_3_0.rad45" "rawButton" vTcl:WidgetProc "Toplevel1" 1
    radiobutton $site_3_0.rad46 \
        -activebackground {#d9d9d9} -activeforeground black \
        -background {#d9d9d9} -command showImageChange -foreground {#000000} \
        -highlightcolor black -justify left -text Filtered -value filtered \
        -variable whichImage 
    vTcl:DefineAlias "$site_3_0.rad46" "filteredButton" vTcl:WidgetProc "Toplevel1" 1
    radiobutton $site_3_0.rad47 \
        -activebackground {#d9d9d9} -activeforeground black \
        -background {#d9d9d9} -command showImageChange -foreground {#000000} \
        -highlightcolor black -justify left -text None -value none \
        -variable whichImage 
    vTcl:DefineAlias "$site_3_0.rad47" "noneButton" vTcl:WidgetProc "Toplevel1" 1
    checkbutton $site_3_0.che48 \
        -activebackground {#d9d9d9} -activeforeground black \
        -background {#d9d9d9} -command invertImageChange \
        -foreground {#000000} -highlightcolor black -justify left \
        -text Invert -variable invertImage 
    vTcl:DefineAlias "$site_3_0.che48" "invertCheck" vTcl:WidgetProc "Toplevel1" 1
    place $site_3_0.rad45 \
        -in $site_3_0 -x 10 -y 20 -anchor nw -bordermode ignore 
    place $site_3_0.rad46 \
        -in $site_3_0 -x 60 -y 20 -anchor nw -bordermode ignore 
    place $site_3_0.rad47 \
        -in $site_3_0 -x 150 -y 20 -anchor nw -bordermode ignore 
    place $site_3_0.che48 \
        -in $site_3_0 -x 210 -y 20 -anchor nw -bordermode ignore 
    labelframe $top.lab41 \
        -foreground black -text Annotations -background {#d9d9d9} -height 45 \
        -highlightcolor black -width 430 
    vTcl:DefineAlias "$top.lab41" "annotationsframe" vTcl:WidgetProc "Toplevel1" 1
    set site_3_0 $top.lab41
    checkbutton $site_3_0.che42 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -command changeVisibleAnnotations \
        -foreground {#000000} -highlightcolor black -justify left \
        -text circles -variable showCircles 
    vTcl:DefineAlias "$site_3_0.che42" "circlesCheck" vTcl:WidgetProc "Toplevel1" 1
    checkbutton $site_3_0.che43 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -command changeVisibleAnnotations \
        -foreground {#000000} -highlightcolor black -justify left \
        -text triang -variable showTriang 
    vTcl:DefineAlias "$site_3_0.che43" "triangulationCheck" vTcl:WidgetProc "Toplevel1" 1
    checkbutton $site_3_0.che49 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -command changeVisibleAnnotations \
        -foreground {#000000} -highlightcolor black -justify left \
        -text defects -variable showDefects 
    vTcl:DefineAlias "$site_3_0.che49" "defectsCheck" vTcl:WidgetProc "Toplevel1" 1
    checkbutton $site_3_0.che50 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -command changeVisibleAnnotations \
        -foreground {#000000} -highlightcolor black -justify left -text angle \
        -variable showOrientation 
    vTcl:DefineAlias "$site_3_0.che50" "orientationCheck" vTcl:WidgetProc "Toplevel1" 1
    checkbutton $site_3_0.che54 \
        -activebackground {#d9d9d9} -activeforeground black -anchor w \
        -background {#d9d9d9} -command changeVisibleAnnotations \
        -foreground {#000000} -highlightcolor black -justify left \
        -text trajectories -variable showTraject 
    vTcl:DefineAlias "$site_3_0.che54" "trajectCheck" vTcl:WidgetProc "Toplevel1" 1
    place $site_3_0.che42 \
        -in $site_3_0 -x 10 -y 20 -anchor nw -bordermode ignore 
    place $site_3_0.che43 \
        -in $site_3_0 -x 90 -y 20 -anchor nw -bordermode ignore 
    place $site_3_0.che49 \
        -in $site_3_0 -x 160 -y 20 -anchor nw -bordermode ignore 
    place $site_3_0.che50 \
        -in $site_3_0 -x 240 -y 20 -anchor nw -bordermode ignore 
    place $site_3_0.che54 \
        -in $site_3_0 -x 310 -y 20 -width 116 -relwidth 0 -height 17 \
        -relheight 0 -anchor nw -bordermode ignore 
    checkbutton $top.che51 \
        -activebackground {#d9d9d9} -activeforeground black \
        -background {#d9d9d9} -command showStatsWin -foreground {#000000} \
        -highlightcolor black -justify left -text {Show stats} \
        -variable showStats -wraplength 40 
    vTcl:DefineAlias "$top.che51" "statsCheck" vTcl:WidgetProc "Toplevel1" 1
    canvas $top.can53 \
        -background {#d9d9d9} -closeenough 1.0 -height 361 \
        -highlightcolor black -insertbackground black -relief ridge \
        -selectbackground {#c4c4c4} -selectforeground black -width 841 
    vTcl:DefineAlias "$top.can53" "imgCanvas" vTcl:WidgetProc "Toplevel1" 1
    bind $top.can53 <B1-Motion> {
        lambda e: xxx(e)
    }
    bind $top.can53 <Button-3> {
        lambda e: xxx(e)
    }
    bind $top.can53 <ButtonRelease-3> {
        lambda e: xxx(e)
    }
    bind $top.can53 <MouseWheel> {
        lambda e: xxx(e)
    }
    ###################
    # SETTING GEOMETRY
    ###################
    place $top.lab40 \
        -in $top -x 10 -y 10 -width 280 -relwidth 0 -height 45 -relheight 0 \
        -anchor nw -bordermode ignore 
    place $top.lab41 \
        -in $top -x 310 -y 10 -width 430 -relwidth 0 -height 45 -relheight 0 \
        -anchor nw -bordermode ignore 
    place $top.che51 \
        -in $top -x 750 -y 20 -anchor nw -bordermode ignore 
    place $top.can53 \
        -in $top -x 10 -y 60 -width 800 -relwidth 0 -height 800 -relheight 0 \
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

