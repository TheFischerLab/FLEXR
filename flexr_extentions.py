import os
import csv
import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GObject, GdkPixbuf
from gi.repository import Gio
from gi.repository import GLib
import coot
import coot_gui_api
import coot_utils
#import coot_gui
from coot_gui import add_simple_action_to_menu, attach_module_menu_button,molecule_chooser_gui_generic,make_store_for_model_molecule_combobox,make_store_for_map_molecule_combobox
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

#---------------------------------------------------------------------
#  FLEXR interface
#---------------------------------------------------------------------

def add_module_flexr():
    menu = attach_module_menu_button('FLEXR')
    add_simple_action_to_menu(
    menu, "FLEXR","run_flexr_gui",lambda _simple_action, _arg: add_flexr_gui())
    add_simple_action_to_menu(
    menu, "FLEX-CHECK","run_multiconfvalidation_gui",lambda _simple_action, _arg: add_multiconfvalidation_gui())
    add_simple_action_to_menu(
    menu, "SWITCHER","coot_switch_alts_button",lambda _simple_action, _arg: add_switch_alts_button())

def run_flexr(imol,mol,imap,mtz,ringerfile,branching,densitythreshold,singleconf,ringerplots,altlimit,densityfilter,clashfilter):

    from flexr import main
    from src.building import flexr_build
    from src.building import flexr_analysis

    print('Initiating FLEXR...')

    from src.flexrpkg.top_level import args

    ARGS = args()


    if (ringerfile == '/') & (imap is not None):
        try:
            print('Running Ringer...')
            print(mol,mtz)

            pdbformat = mol[-4:]

            # use phenix to remove alt confs
            os.system('phenix.pdbtools %s remove_alt_confs=True' % (mol))

            if pdbformat == '.ent':
                pdbformat = '.pdb'
                mol = mol[:-4]+pdbformat

            mol = mol.split('/')[-1][:-4]+"_modified"+pdbformat
            imolnoconf = coot.read_pdb(mol)
            ARGS.cootmolnum = imolnoconf

            os.system('phenix.maps %s %s' % (mol,mtz))
            map_coef = mol.split()[0][:-4]+"_map_coeffs.mtz"
            #print(map_coef)

            os.system('mmtbx.ringer %s %s sampling_angle=2' % (mol,map_coef))
            ringerfile = mol.split('/')[-1][:-4]+"_ringer.csv"
            ARGS.filename = ringerfile
            altsfile = ringerfile[:-4]+"_"+str(densitythreshold)+"_alts.csv"

            testfileexists = open(ringerfile)
            testfileexists.close()

        except:
            print('Cannot find Phenix/CCTBX or it failed.')
            print('Done.')
            raise ValueError

    else:
        ARGS.filename = ringerfile
        altsfile = ringerfile[:-4]+"_"+str(densitythreshold)+'_alts.csv'

    # options in the gui
    ARGS.sigmathreshold = densitythreshold
    ARGS.build_limit = altlimit
    ARGS.branching = branching
    ARGS.singleconfs = singleconf
    ARGS.plot = ringerplots
    ARGS.step = 2
    ARGS.clashscore = clashfilter
    ARGS.densityscore = densityfilter

    # these opts so we can call the build function instead of running the script.
    ARGS.pdb = 'None'
    ARGS.build = False

    # run the main flexr script
    main(ARGS)

    # run the building step
    flexrmolnum = flexr_build.building_run(altsfile,ARGS.pdb,ARGS.branching,ARGS.cootmolnum,ARGS.exitcoot,ARGS.clashscore,ARGS.densityscore)
    #try:
    flexr_analysis.output_summaries(altsfile,imol,flexrmolnum)
    #except:
    #    print('Error producing output summary.')

    print('FLEXR is finished.')
    print('')

def add_flexr_gui():

    def delete_event(*args):
        window.destroy()
        return False

    # get current MODEL for building
    def get_molecule():
        tree_iter = combobox_molecule.get_active_iter()
        imol = -1
        if tree_iter is not None:
            model = combobox_molecule.get_model()
            it = model[tree_iter]
            imol = it[0]
            name = coot.molecule_name(imol)
        return imol,name

    # get current MAP for building
    def get_mtz():
        tree_iter = combobox_map.get_active_iter()
        imap = -1
        if tree_iter is not None:
            mtz = combobox_map.get_model()
            it = mtz[tree_iter]
            imap = it[0]
            name = coot.molecule_name(imap)
        return imap,name

    # some user input, text
    def get_ringer_path():
        t = dist_entry.get_text()
        try:
            ret = str(t)
        except:
            ret = False
        return ret

    # some user input, text
    def get_threshold():
        t = dist_entry2.get_text()
        try:
            ret = float(t)
        except:
            ret = False
        return ret

    # some user input, text
    alt_limit_list = list(range(2,11))
    def get_alt_limit():
        t = alt_limit.get_active()
        try:
            ret = int(t)
            ret = alt_limit_list[ret]
        except:
            ret = False
        return ret

    # Sidechain branching
    # some user input, drop down box
    def get_getbranching():
        at = combobox_coordination.get_active()
        n = str(at)
        return n

    # Single conformer modelling
    # some user input, drop down box
    def get_singleconf():
        at = combobox_coordination2.get_active()
        n = bool(at)
        return n

    # Ringer plots
    # some user input, drop down box
    def get_ringerplots():
        at = combobox_coordination3.get_active()
        n = bool(at)
        return n

    #filtering options
    def get_filter_ops():
        if clashcheck_button.get_active():
            clashfilter = True
        else:
            clashfilter = False
        if densitycheck_button.get_active():
            densityfilter = True
        else:
            densityfilter = False
        return clashfilter,densityfilter

    def apply_cb(*args):
        imol,mol = get_molecule()
        imap,mtz = get_mtz()
        ringerpath = get_ringer_path()
        densitythreshold = get_threshold()
        altlimit = get_alt_limit()
        branching = get_getbranching()
        singleconf = get_singleconf()
        ringerplots = get_ringerplots()
        clashfilter,densityfilter = get_filter_ops()
        if ringerpath:
            run_flexr(imol,mol,imap,mtz,ringerpath,branching,densitythreshold,singleconf,ringerplots,altlimit,densityfilter,clashfilter)
            #print("now call update_water_results() with imol", imol, "coordination number", n, "imol", imol)
            #update_water_results(imol, n, d)

    # design GUI

    window = Gtk.Window()
    window.set_title("FLEXR: multi-conformer modeling")
    window.set_default_size(300,200)
    vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
    h_sep = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)

    frame = Gtk.Frame()
    frame.set_size_request(100, 200)

    #header image
    from src.flexrpkg.top_level import get_coot_loc
    libraryloc, cootloc,cootexe = get_coot_loc()
    imagepath=cootloc+"/img/logo.png"
    picture = Gtk.Picture.new_for_filename(imagepath)
    picture.set_content_fit(Gtk.ContentFit.CONTAIN)
    #picture.set_can_shrink(True)
    picture.set_size_request(100, 200)

    frame.set_child(picture)
    vbox.append(frame)

    ## header
    hint_text = Gtk.Label()
    hint_text.set_markup("\n For full instructions, go to our <a href=\"https://github.com/TheFischerLab/FLEXR-GUI\" "
                 "title=\"\">GitHub</a> repo")
    vbox.append(hint_text)

    hint_text2 = Gtk.Label()
    hint_text2.set_markup("\n Please cite: <a href=\"https://doi.org/10.1107/S2059798323002498\" "
                "title=\"\">Stachowski and Fischer, Acta Cryst. (2023) D79, 354-357</a> \n")
    vbox.append(hint_text2)

    ## drop down for picking molecule
    hint_text = Gtk.Label(label="\n Molecule:")
    combobox_molecule = Gtk.ComboBox()
    vbox.append(hint_text)
    combobox_mol_items = make_store_for_model_molecule_combobox()
    combobox_molecule.set_model(combobox_mol_items)

    renderer_text = Gtk.CellRendererText()
    if len(combobox_mol_items) > 0:
        combobox_molecule.set_active(0)
    combobox_molecule.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
    combobox_molecule.pack_start(renderer_text, True)
    combobox_molecule.add_attribute(renderer_text, "text", 1)
    combobox_molecule.set_margin_start(5)
    combobox_molecule.set_margin_end(5)
    combobox_molecule.set_margin_top(10)
    combobox_molecule.set_margin_bottom(10)
    vbox.append(combobox_molecule)

    # some debugging stuff
    print("debug:: add_flexr_gui(): combobox_mol_items:", combobox_mol_items)
    print("debug:: add_flexr_gui(): combobox_molecule:",  combobox_molecule)

    combobox_molecule.set_active(0)

    ## drop down for picking map
    hint_text = Gtk.Label(label="\n MTZ:")
    combobox_map = Gtk.ComboBox()
    vbox.append(hint_text)
    combobox_mol_items2 = make_store_for_map_molecule_combobox()
    combobox_map.set_model(combobox_mol_items2)

    renderer_text = Gtk.CellRendererText()
    if len(combobox_mol_items2) > 0:
        combobox_map.set_active(0)
    combobox_map.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
    combobox_map.pack_start(renderer_text, True)
    combobox_map.add_attribute(renderer_text, "text", 1)
    combobox_map.set_margin_start(5)
    combobox_map.set_margin_end(5)
    combobox_map.set_margin_top(10)
    combobox_map.set_margin_bottom(10)
    vbox.append(combobox_map)

    # some debugging stuff
    print("debug:: add_flexr_gui(): combobox_mol_items:", combobox_mol_items)
    print("debug:: add_flexr_gui(): combobox_molecule:",  combobox_molecule)

    combobox_map.set_active(0)


    results_vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
    def clear_previous_results():
        for this_vbox in [results_vbox, metal_results_vbox]:
            child = this_vbox.get_first_child()
            while child is not None:
                next_child = child.get_next_sibling()
                this_vbox.remove(child)
                child = next_child

    window.set_child(vbox)

    # default values for text entry
    hbox_max_dist = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    #dist_label = Gtk.Label(label="FLEXR alt file location: ")
    dist_label = Gtk.Label(label="Ringer CSV (if precomputed): ")
    dist_entry = Gtk.Entry()
    hbox_max_dist.append(dist_label)
    hbox_max_dist.append(dist_entry)
    hbox_max_dist.set_margin_start(6)
    hbox_max_dist.set_margin_end(6)
    hbox_max_dist.set_margin_top(4)
    hbox_max_dist.set_margin_bottom(4)
    vbox.append(hbox_max_dist)
    dist_entry.set_text("/")
    #dist_entry.set_width_chars(5)

    # default values for text entry
    density_threshold = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    density_threshold_label = Gtk.Label(label="Density threshold: ")
    dist_entry2 = Gtk.Entry()
    density_threshold.append(density_threshold_label)
    density_threshold.append(dist_entry2)
    density_threshold.set_margin_start(6)
    density_threshold.set_margin_end(6)
    density_threshold.set_margin_top(4)
    density_threshold.set_margin_bottom(4)
    vbox.append(density_threshold)
    dist_entry2.set_text("0.30")
    #dist_entry.set_width_chars(5)

    # default values for text entry
    #alt_limit = Gtk.ComboBoxText(orientation=Gtk.Orientation.HORIZONTAL)
    alt_limit = Gtk.ComboBoxText.new()
    alt_limit_label = Gtk.Label(label="  Max. number of alts/residue: ")
    #dist_entry3 = Gtk.Entry()
    #alt_limit.append(alt_limit_label)
    #alt_limit.append(dist_entry3)
    for i in alt_limit_list:
        alt_limit.append_text(str(i))
    alt_limit.set_active(1)
    alt_limit.set_margin_start(6)
    alt_limit.set_margin_end(1)
    alt_limit.set_margin_top(4)
    alt_limit.set_margin_bottom(4)
    hbox_chooser = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    hbox_number_chooser = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    vbox.append(hbox_chooser)
    hbox_number_chooser.append(alt_limit_label)
    hbox_number_chooser.append(alt_limit)
    vbox.append(hbox_number_chooser)
    #dist_entry3.set_text("3")
    #dist_entry.set_width_chars(5)

    # Create places for drop down option
    # coordination number combobox
    number_text = Gtk.Label(label="  Sidechain branching: ")
    combobox_coordination = Gtk.ComboBoxText.new()
    for i in ['ALL','CA']:
        combobox_coordination.append_text(str(i))

    combobox_coordination.set_active(0)
    combobox_coordination.set_margin_start(6)
    combobox_coordination.set_margin_end(6)
    combobox_coordination.set_margin_top(4)
    combobox_coordination.set_margin_bottom(4)
    hbox_chooser = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    hbox_number_chooser = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    vbox.append(hbox_chooser)
    hbox_number_chooser.append(number_text)
    hbox_number_chooser.append(combobox_coordination)
    vbox.append(hbox_number_chooser)

    # Create places for drop down option
    # coordination number combobox
    number_text2 = Gtk.Label(label="  Single conformer modelling (beta): ")
    combobox_coordination2 = Gtk.ComboBoxText.new()
    for i in ['False','True']:
        combobox_coordination2.append_text(str(i))

    combobox_coordination2.set_active(0)
    combobox_coordination2.set_margin_start(6)
    combobox_coordination2.set_margin_end(6)
    combobox_coordination2.set_margin_top(4)
    combobox_coordination2.set_margin_bottom(4)
    hbox_chooser2 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    hbox_number_chooser2 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    vbox.append(hbox_chooser2)
    hbox_number_chooser2.append(number_text2)
    hbox_number_chooser2.append(combobox_coordination2)
    vbox.append(hbox_number_chooser2)

    # Create places for drop down option
    # coordination number combobox
    number_text3 = Gtk.Label(label="  Ringer plotting (slow): ")
    combobox_coordination3 = Gtk.ComboBoxText.new()
    for i in ['False','True']:
        combobox_coordination3.append_text(str(i))

    combobox_coordination3.set_active(0)
    combobox_coordination3.set_margin_start(6)
    combobox_coordination3.set_margin_end(6)
    combobox_coordination3.set_margin_top(4)
    combobox_coordination3.set_margin_bottom(4)
    hbox_chooser3 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    hbox_number_chooser3 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    vbox.append(hbox_chooser3)
    hbox_number_chooser3.append(number_text3)
    hbox_number_chooser3.append(combobox_coordination3)
    vbox.append(hbox_number_chooser3)

    # Create places for drop down option
    # coordination number combobox
    number_text4 = Gtk.Label(label="Filter alts by clash?        ")
    hbox_chooser4 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    hbox_chooser4.set_margin_start(6)
    hbox_chooser4.set_margin_bottom(4)
    hbox_chooser4.set_margin_top(4)
    clashcheck_button_local = Gtk.CheckButton(label = "")
    clashcheck_button = clashcheck_button_local
    clashcheck_button.set_active(True)

    vbox.append(hbox_chooser4)
    hbox_chooser4.append(number_text4)
    hbox_chooser4.append(clashcheck_button_local)


    #
    number_text5 = Gtk.Label(label="Filter alts by density?     ")
    hbox_chooser5 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    hbox_chooser5.set_margin_start(6)
    hbox_chooser5.set_margin_bottom(4)
    hbox_chooser5.set_margin_top(4)

    densitycheck_button_local = Gtk.CheckButton(label = "")
    densitycheck_button = densitycheck_button_local
    densitycheck_button.set_active(False)

    vbox.append(hbox_chooser5)
    hbox_chooser5.append(number_text5)
    hbox_chooser5.append(densitycheck_button_local)

    ## default values for text entry
    #density_threshold = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    #density_threshold_label = Gtk.Label(label="Density threshold: ")
    #dist_entry2 = Gtk.Entry()
    #density_threshold.append(density_threshold_label)
    #density_threshold.append(dist_entry2)
    #density_threshold.set_margin_start(6)
    #density_threshold.set_margin_end(6)
    #density_threshold.set_margin_top(4)
    #density_threshold.set_margin_bottom(4)
    #vbox.append(density_threshold)
    #dist_entry2.set_text("0.30")
    ##dist_entry.set_width_chars(5)


    #combobox_coordination4 = Gtk.ComboBoxText.new()
    #radiobutton_yes = Gtk.ToggleButton(None,True)
    #radiobutton_no = Gtk.ToggleButton(radiobutton_yes,False)
    #radiobutton_yes.set_active(True)
    #radiobutton_yes.show()
    #radiobutton_no.show()
    #combobox_coordination4.append(radiobutton_yes)
    #combobox_coordination4.append(radiobutton_no)
#
    #combobox_coordination4.set_active(0)
    #combobox_coordination4.set_margin_start(6)
    #combobox_coordination4.set_margin_end(6)
    #combobox_coordination4.set_margin_top(4)
    #combobox_coordination4.set_margin_bottom(4)
    #hbox_chooser4 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    #hbox_number_chooser4 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    #vbox.append(hbox_chooser4)
    #hbox_number_chooser4.append(number_text4)
    #hbox_number_chooser4.append(combobox_coordination4)
    #vbox.append(hbox_number_chooser3)

    ## Place for 'Output 2'
    #scrolled_win = Gtk.ScrolledWindow()
    #metal_results_scrolled_win = Gtk.ScrolledWindow()
    #metal_results_frame = Gtk.Frame()
    #metal_results_vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
    #metal_results_label = Gtk.Label(label="Output 1")
    #metal_results_scrolled_win.set_child(metal_results_frame)
    #metal_results_frame.set_child(metal_results_vbox)
    #vbox.append(metal_results_label)
    #vbox.append(metal_results_scrolled_win)

    ## Place for 'Output 1'
    #water_results_label = Gtk.Label(label="Output 2")
    #scrolled_win.set_child(results_vbox)
    #vbox.append(water_results_label)
    #vbox.append(scrolled_win)
    #vbox.append(h_sep)

    ## Place exit/run buttons
    apply_button = Gtk.Button(label="  Run  ")
    close_button = Gtk.Button(label="  Close  ")
    hbox_buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)

    # create buttons to exit/run
    apply_button.set_margin_start(10)
    apply_button.set_margin_end(10)
    apply_button.set_margin_top(4)
    apply_button.set_margin_bottom(4)
    close_button.set_margin_start(6)
    close_button.set_margin_end(6)
    close_button.set_margin_top(4)
    close_button.set_margin_bottom(4)
    hbox_buttons.set_margin_start(6)
    hbox_buttons.set_margin_end(6)
    hbox_buttons.set_margin_top(4)
    hbox_buttons.set_margin_bottom(4)
    hbox_buttons.append(close_button)
    hbox_buttons.append(apply_button)
    vbox.append(hbox_buttons)
    close_button.connect("clicked", delete_event)
    #dist_entry.connect("activate", entry_activate_event)
    apply_button.connect("clicked", apply_cb)

    window.show()

#---------------------------------------------------------------------
#  FLEX-CHECK interface
#---------------------------------------------------------------------

def run_validation(imol,mol):

    ## convert cif to pdb
    pdbformat = mol[-4:]
    if pdbformat == ".cif":
        coot.write_pdb_file(imol,mol[:-4]+"_flex_check.pdb")
        mol = mol[:-4]+"_flex_check.pdb"
    else:
        mol = mol[:-4]+pdbformat

    print('')
    print('')
    print('Starting validation from Coot...')
    print('')
    print('Current working directory: ')
    os.system('pwd')
    print('Input: %s' % (mol))
    print('')
    print('')

    ## try to find installation
    from src.flexrpkg.top_level import get_coot_loc
    libraryloc, cootloc,cootexe = get_coot_loc()
    try:
        buildpath = cootloc+'/src/tools/flex_check.py'
        os.system('/Applications/ccp4-9/coot_py3/Frameworks/Python.framework/Versions/3.9/bin/python3.9 %s %s' % (buildpath,mol))
    except:
        print('Cannot find FLEXR installation.')

    if os.path.exists('multiconf_refinement_check_output.csv'):
        os.replace('multiconf_refinement_check_output.csv','./validation/multiconf_refinement_check_output.csv')

def add_multiconfvalidation_gui():

    def delete_event(*args):
        window.destroy()
        return False

    # get current MODEL for building
    def get_molecule():
        tree_iter = combobox_molecule.get_active_iter()
        imol = -1
        if tree_iter is not None:
            model = combobox_molecule.get_model()
            it = model[tree_iter]
            imol = it[0]
            name = coot.molecule_name(imol)
            #coot.turn_off_backup(imol)
        return imol,name

    def apply_cb(*args):
        imol,mol = get_molecule()
        run_validation(imol,mol)

        #update CSV view
        last = vbox.get_last_child()
        vbox.remove(last)
        treeview = Gtk.TreeView()
        treeview.set_grid_lines(Gtk.TreeViewGridLines.BOTH)
        treeview.set_vexpand(True)
        scrolled = Gtk.ScrolledWindow()
        scrolled.set_child(treeview)
        scrolled.set_vexpand(True)
        load_csv_into_treeview(treeview, './validation/multiconf_refinement_check_output.csv')
        vbox.append(scrolled)
        #

        treeview.connect("row-activated", coot_view_csv_row)

    def remove_lonely_Hs(*args):
        imol,mol = get_molecule()
        with open('./validation/multiconf_refinement_check_output.csv', newline='') as csvfile:
            reader = csv.reader(csvfile)
            rows = [row for row in reader]
        alts = ['A','B','C','D','E','F']
        for row in rows[1:]:
            if 'X' in row[8]:
                for alt in alts:
                    print('baddie!')
                    coot.delete_atom(imol,row[1],int(row[2]),'',' H  ',alt)

        print("Lonely Hydrogens Removed.")
        print("Done.")

    def match_occupancies(*args):
        imol,mol = get_molecule()
        with open('./validation/multiconf_refinement_check_output.csv', newline='') as csvfile:
            reader = csv.reader(csvfile)
            rows = [row for row in reader]
        for row in rows[1:]:
            if ('X' in row[7]) or ('X' in row[9]):
                res_info = coot.residue_info_py(imol,row[1],int(row[2]),'')
                alts = []
                for atom in res_info:
                    alt = atom[0][1]
                    alts.append(alt)
                alts = list(set(alts))
                if (len(alts)>1):
                    if ('' in alts):
                        new_occ = 1/(len(alts)-1)
                    else:
                        new_occ = 1/(len(alts))
                #new_occ = "{:.2f}".format(new_occ)
                for j in res_info:
                    atom = j[0][0]
                    alt = j[0][1]
                    if alt == '':
                        #update occ
                        coot.set_atom_attribute(imol,row[1],int(row[2]),'',atom,alt,'occ',float(1.0))
                    else:
                        coot.set_atom_attribute(imol,row[1],int(row[2]),'',atom,alt,'occ',new_occ)
        print("Mismatched occupancies reset.")
        print("Done.")


    # design GUI

    window = Gtk.Window()
    window.set_title("FLEX-check: multi-conformer model validation")
    window.set_default_size(300,625)
    vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
    h_sep = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)

    frame = Gtk.Frame()
    frame.set_size_request(100, 125)

    #header image
    from src.flexrpkg.top_level import get_coot_loc
    libraryloc, cootloc,cootexe = get_coot_loc()
    imagepath=cootloc+"/img/check.png"
    picture = Gtk.Picture.new_for_filename(imagepath)
    #picture.set_content_fit(Gtk.ContentFit.CONTAIN)
    #picture.set_can_shrink(True)
    picture.set_size_request(100, 125)

    frame.set_child(picture)
    vbox.append(frame)

    ## header
    hint_text = Gtk.Label()
    hint_text.set_markup("\n For full instructions, go to our <a href=\"https://github.com/TheFischerLab/FLEXR-GUI\" "
                 "title=\"\">GitHub</a> repo")
    vbox.append(hint_text)

    hint_text2 = Gtk.Label()
    hint_text2.set_markup("\n Please cite: <a href=\"https://doi.org/10.1107/S2059798323002498\" "
                "title=\"\">Stachowski and Fischer, Acta Cryst. (2023) D79, 354-357</a> \n")
    vbox.append(hint_text2)

    ## drop down for picking molecule
    hint_text = Gtk.Label(label="\n Molecule:")
    combobox_molecule = Gtk.ComboBox()
    vbox.append(hint_text)
    combobox_mol_items = make_store_for_model_molecule_combobox()
    combobox_molecule.set_model(combobox_mol_items)

    renderer_text = Gtk.CellRendererText()
    if len(combobox_mol_items) > 0:
        combobox_molecule.set_active(0)
    combobox_molecule.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
    combobox_molecule.pack_start(renderer_text, True)
    combobox_molecule.add_attribute(renderer_text, "text", 1)
    combobox_molecule.set_margin_start(5)
    combobox_molecule.set_margin_end(5)
    combobox_molecule.set_margin_top(10)
    combobox_molecule.set_margin_bottom(10)
    vbox.append(combobox_molecule)

    # some debugging stuff
    print("debug:: add_multiconfvalidation_gui(): combobox_mol_items:", combobox_mol_items)
    print("debug:: add_multiconfvalidation_gui(): combobox_molecule:",  combobox_molecule)

    combobox_molecule.set_active(0)


    results_vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
    def clear_previous_results():
        for this_vbox in [results_vbox, metal_results_vbox]:
            child = this_vbox.get_first_child()
            while child is not None:
                next_child = child.get_next_sibling()
                this_vbox.remove(child)
                child = next_child

    window.set_child(vbox)

    ## Place exit/run buttons
    apply_button = Gtk.Button(label="  Run  ")
    close_button = Gtk.Button(label="  Close  ")
    clean_h_button = Gtk.Button(label="  Clean Hs  ")
    match_occs_button = Gtk.Button(label="  Reset Occ  ")
    switch_alts_button = Gtk.Button(label="  Switch Alts  ")

    hbox_buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    hbox_buttons.set_halign(Gtk.Align.CENTER)

    hbox_buttons2 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    hbox_buttons2.set_halign(Gtk.Align.CENTER)

    # create buttons to exit/run
    apply_button.set_margin_start(6)
    apply_button.set_margin_end(6)
    apply_button.set_margin_top(4)
    apply_button.set_margin_bottom(4)

    close_button.set_margin_start(6)
    close_button.set_margin_end(6)
    close_button.set_margin_top(4)
    close_button.set_margin_bottom(4)

    clean_h_button.set_margin_start(2)
    clean_h_button.set_margin_end(2)
    clean_h_button.set_margin_top(4)
    clean_h_button.set_margin_bottom(4)

    match_occs_button.set_margin_start(2)
    match_occs_button.set_margin_end(2)
    match_occs_button.set_margin_top(4)
    match_occs_button.set_margin_bottom(4)

    switch_alts_button.set_margin_start(2)
    switch_alts_button.set_margin_end(2)
    switch_alts_button.set_margin_top(4)
    switch_alts_button.set_margin_bottom(4)

    #hbox_buttons.set_margin_start(6)
    #hbox_buttons.set_margin_end(6)
    #hbox_buttons.set_margin_top(4)
    #hbox_buttons.set_margin_bottom(4)

    hbox_buttons.append(close_button)
    hbox_buttons.append(apply_button)


    hbox_buttons2.append(clean_h_button)
    hbox_buttons2.append(match_occs_button)
    hbox_buttons2.append(switch_alts_button)

    vbox.append(hbox_buttons)
    vbox.append(hbox_buttons2)
    close_button.connect("clicked", delete_event)
    #dist_entry.connect("activate", entry_activate_event)


    # viewer

    # CSV Table
    treeview = Gtk.TreeView()
    treeview.set_grid_lines(Gtk.TreeViewGridLines.BOTH)
    treeview.set_vexpand(True)
    scrolled = Gtk.ScrolledWindow()
    scrolled.set_child(treeview)
    scrolled.set_vexpand(True)


    if not os.path.exists('./validation/multiconf_refinement_check_output.csv'):
        #open('./validation/multiconf_refinement_check_output.csv', 'a').close()
        os.mkdir('validation')
        with open('./validation/multiconf_refinement_check_output.csv','w+') as f:
            f.write('model,chain,res_num,res_type,alt_loc,occupancy,1. Occ<0.1,2. MismatchedOcc,3. Lonely_Hs,4. SumOcc!=1,5. Wrong Highest Alt ID,6. Wrong # of Alt IDs')
        f.close()

#    if os.path.exists('./validation/multiconf_refinement_check_output.csv'):
#        with open('./validation/multiconf_refinement_check_output.csv', newline='') as csvfile:
#            load_csv_into_treeview(treeview, './validation/multiconf_refinement_check_output.csv')

    if os.path.exists('./validation/multiconf_refinement_check_output.csv'):
        #open('./validation/multiconf_refinement_check_output.csv', 'a').close()
        #os.mkdir('validation')
        with open('./validation/multiconf_refinement_check_output.csv','w+') as f:
            f.write('model,chain,res_num,res_type,alt_loc,occupancy,1. Occ<0.1,2. MismatchedOcc,3. Lonely_Hs,4. SumOcc!=1,5. Wrong Highest Alt ID,6. Wrong # of Alt IDs')
        f.close()

    vbox.append(scrolled)
    apply_button.connect("clicked", apply_cb,vbox)

    clean_h_button.connect("clicked", remove_lonely_Hs)
    match_occs_button.connect("clicked", match_occupancies)
    switch_alts_button.connect("clicked", add_switch_alts_button)


    window.show()

def load_csv_into_treeview(treeview, csv_path):
    # Remove previous columns
    for col in treeview.get_columns():
        treeview.remove_column(col)
    # Load CSV
    with open(csv_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        rows = [row for row in reader]
    if not rows:
        return
    num_columns = len(rows[0])
    types = [str] * num_columns
    liststore = Gtk.ListStore(*types)
    for row in rows[1:]:
        liststore.append(row)
    treeview.set_model(liststore)
    for i, col_title in enumerate(rows[0]):
        renderer = Gtk.CellRendererText()
        column = Gtk.TreeViewColumn(col_title, renderer, text=i)
        column.set_cell_data_func(renderer, cell_background_func)
        treeview.append_column(column)
    return rows[0]

def coot_view_csv_row(treeview,path,column):
    model = treeview.get_model()
    iter = model.get_iter(path)
    row_data = model[iter][:]
    print("Row clicked:", row_data)
    coot.set_go_to_atom_chain_residue_atom_name(row_data[1],int(row_data[2]),'CA')

def cell_background_func(column, cell, model, iter, data):
    path = model.get_path(iter)
    row_index = path.get_indices()[0]
    color = "#a9a9a9" if row_index % 2 == 0 else "dimgray"
    cell.set_property("cell-background", color)

#---------------------------------------------------------------------
#  Coot Switch Alt IDs interface
#---------------------------------------------------------------------

def add_module_switch_alts_button():
    menu = attach_module_menu_button('SWITCHER')
    add_simple_action_to_menu(
    menu, "Open","coot_switch_alts_button",lambda _simple_action, _arg: add_switch_alts_button())

def add_switch_alts_button(*args):

    def delete_event(*args):
        window.destroy()
        return False

    def get_event(*args):
        try:
            #get info from click
            atom_info = coot.select_atom_under_pointer_py()
            imol = atom_info[0]
            res_info = atom_info[1]
            chain = res_info[1]
            resn = res_info[2]
            ins_code = res_info[3]
            atom = res_info[4]
            alt = res_info[5]

            res_info = coot.residue_info_py(imol,chain,resn,ins_code)
            #print(res_info)
            cas = []
            for i in range(len(res_info)):
                atom = res_info[i][0]
                info = res_info[i][1]
                coords = res_info[i][2]
                #if atom[0] == ' CB ':
                if atom[1] != '':
                    cas.append([str(info[0]),atom[1]])
            cas = list(set(map(tuple,cas)))
            if len(cas) == 2:
                cas.append(['',''])
            #print(cas)
            #return imol,res_info,chain,resn,ins_code,atom_alt
            #print(len(entries))
            if len(cas) == 2:
                cas.append(['',''])
            if len(cas) > 0:
                casT = [list(row) for row in zip(*cas)]
                for row in [0,1]:
                    for column in [0,1,2]:
                        entries[row][column].set_text(casT[row][column])

            return res_info,cas

        except:
            print('No alts at this residue')
            for row in [0,1]:
                for column in [0,1,2]:
                    entries[row][column].set_text('')

    def update_event(*args):

        new_data = []
        for row in [0,1]:
            for column in [0,1,2]:
                t = entries[row][column].get_text()
                new_data.append(t)

        tuple_data = [[new_data[0],new_data[3]],[new_data[1],new_data[4]],[new_data[2],new_data[5]]]

        ## shouldn't have to re calculate...
        #get info from click
        atom_info = coot.select_atom_under_pointer_py()
        imol = atom_info[0]
        res_info = atom_info[1]
        chain = res_info[1]
        resn = res_info[2]
        ins_code = res_info[3]
        atom = res_info[4]
        alt = res_info[5]

        res_info = coot.residue_info_py(imol,chain,resn,ins_code)

        cas = []
        for i in range(len(res_info)):
            atom = res_info[i][0]
            info = res_info[i][1]
            coords = res_info[i][2]
            #if atom[0] == ' CB ':
            if atom[1] != '':
                cas.append([str(info[0]),atom[1]])
        cas = list(set(map(tuple,cas)))
        if len(cas) == 2:
            cas.append(['',''])

        for i in [0,1,2]:
            for j in res_info:
                atom = j[0][0]
                alt = j[0][1]
                if (alt == cas[i][1]) & (alt != ''):
                    #update occ
                    coot.set_atom_attribute(imol,chain,resn,ins_code,atom,cas[i][1],'occ',float(tuple_data[i][0]))
                    #update alt
                    coot.set_atom_string_attribute(imol,chain,resn,ins_code,atom,cas[i][1],'alt-conf',tuple_data[i][1]+"z")

        res_info_update = coot.residue_info_py(imol,chain,resn,ins_code)
        for k in res_info_update:
            atom = k[0][0]
            alt = k[0][1]
            #print(atom,alt)
            #print(alt.strip('z'))
            coot.set_atom_string_attribute(imol,chain,resn,ins_code,atom,alt,'alt-conf',alt.strip('z'))

    ### DESIGN GUI ###

    window = Gtk.Window()
    window.set_title("SWITCHER: switch Alt IDs")
    window.set_default_size(300,200)
    vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
    h_sep = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)

    window.set_child(vbox)

    frame = Gtk.Frame()
    frame.set_size_request(110, 103)

    from src.flexrpkg.top_level import get_coot_loc
    libraryloc, cootloc,cootexe = get_coot_loc()
    imagepath=cootloc+"/img/switch.png"

    picture = Gtk.Picture.new_for_filename(imagepath)
    frame.set_child(picture)
    #picture.set_content_fit(Gtk.ContentFit.COVER)
    vbox.append(frame)

    #header
    hint_text = Gtk.Label()
    hint_text.set_markup("1. Double click a residue -> 'Get'\n2. Make changes\n3. Click 'Update'")
    hint_text.set_margin_top(20)
    hint_text.set_margin_bottom(20)
    vbox.append(hint_text)



    # Create a grid
    grid = Gtk.Grid()
    grid.set_row_spacing(10)
    grid.set_column_spacing(10)
    grid.set_margin_top(10)
    grid.set_margin_bottom(20)
    grid.set_margin_start(20)
    grid.set_margin_end(20)

    # Create and attach 4 entry boxes in a 2x2 grid
    row_labels=["Occupancy","Alt ID"]
    entries = []
    for row in range(2):
        row_entries = []
        for col in range(4):
            if col == 0:
                label = Gtk.Label(label=row_labels[row])
                label.set_halign(Gtk.Align.START)
                grid.attach(label,0,row,1,1)
            else:
                entry = Gtk.Entry()
                entry.set_width_chars(1)
                entry.set_max_length(4)
                grid.attach(entry, col, row, 1, 1)
                row_entries.append(entry)
        entries.append(row_entries)

    vbox.append(grid)

    ## Place exit/run buttons
    get_button = Gtk.Button(label="Get")
    update_button = Gtk.Button(label="Update")
    close_button = Gtk.Button(label="Close")

    hbox_buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)

    # create buttons to exit/run
    update_button.set_margin_start(6)
    update_button.set_margin_end(6)
    update_button.set_margin_top(4)
    update_button.set_margin_bottom(4)

    get_button.set_margin_start(6)
    get_button.set_margin_end(6)
    get_button.set_margin_top(4)
    get_button.set_margin_bottom(4)

    close_button.set_margin_start(6)
    close_button.set_margin_end(6)
    close_button.set_margin_top(4)
    close_button.set_margin_bottom(4)

    hbox_buttons.set_margin_start(6)
    hbox_buttons.set_margin_end(6)
    hbox_buttons.set_margin_top(4)
    hbox_buttons.set_margin_bottom(4)

    hbox_buttons.append(close_button)
    hbox_buttons.append(get_button)
    hbox_buttons.append(update_button)

    get_button.connect("clicked", get_event)
    update_button.connect("clicked", update_event)
    close_button.connect("clicked", delete_event)

    vbox.append(hbox_buttons)

    window.show()

#### FINAL ####

add_module_flexr()
