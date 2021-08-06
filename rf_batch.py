#! /usr/bin/env python3
'''
Create a Reflow run file that batch-executes many other Reflow run files.
'''
from os import environ, makedirs, remove
from os.path import abspath, expanduser, isfile
from subprocess import Popen, PIPE
from sys import platform, stderr, stdout
import argparse
import sys

# helpful constants
GUI = False
REFLOW_EXEC_LINUX_URL = 'https://github.com/grailbio/reflow/releases/download/reflow1.6.0/reflow1.6.0.linux.amd64'
REFLOW_EXEC_MAC_URL = 'https://github.com/grailbio/reflow/releases/download/reflow1.6.0/reflow1.6.0.darwin.amd64'

# help text
HELP_TEXT_ABSPATH = "Use absolute paths to Reflow run files"
HELP_TEXT_RUN = "Run batch RF file using 'reflow run'"

# parse user args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output Batch Reflow Run File (RF)")
    parser.add_argument('-a', '--abs_path', action="store_true", help=HELP_TEXT_ABSPATH)
    parser.add_argument('-r', '--run', action="store_true", help=HELP_TEXT_RUN)
    parser.add_argument('rf_files', metavar='RF', type=str, nargs='+', help="Input Reflow Run Files (RF)")
    args = parser.parse_args()
    for rf in args.rf_files:
        if not isfile(rf):
            stderr.write("ERROR: File not found: %s\n" % args.output); exit(1)
    if args.output == 'stdout':
        args.output = stdout
    else:
        if isfile(args.output):
            stderr.write("ERROR: Output file exists: %s\n" % args.output); exit(1)
        args.output = open(args.output, 'w')
    return args

# main function
def main():
    args = parse_args()
    args.output.write("val Main = [\n")
    for rf in args.rf_files:
        if args.abs_path:
            p = abspath(rf)
        else:
            p = rf
        args.output.write('    make("%s").Main,\n' % p)
    args.output.write("]\n")
    args.output.close()
    return args

# run GUI
def run_gui():
    try:
        # imports
        from tkinter import Button, Checkbutton, END, Entry, Frame, IntVar, Label, OptionMenu, StringVar, Tk
        from tkinter.filedialog import askopenfilenames, asksaveasfilename
        from tkinter.scrolledtext import ScrolledText

        # helper function to make a popup
        def gui_popup(message, title=None):
            popup = Tk()
            if title:
                popup.wm_title(title)
            label = Label(popup, text=message)
            label.pack()
            button_close = Button(popup, text="Close", command=popup.destroy)
            button_close.pack(padx=3, pady=3)
            popup.mainloop()

        # create applet
        root = Tk()
        root.geometry("600x450")

        # set up main frame
        frame = Frame(root)
        frame.pack()

        # add header
        header = Label(frame, text="Reflow Batch", font=('Arial',24))
        header.pack()

        # handle output batch RF file selection
        button_out_prefix = "Output Batch RF:\n"
        button_out_nofile = "<none selected>"
        def find_filename_out():
            fn = asksaveasfilename(title="Select Output Batch RF File", filetypes=(('RF file','*.rf'),))
            if len(fn) == 0:
                button_out.configure(text="%s%s" % (button_out_prefix,button_out_nofile))
            else:
                button_out.configure(text="%s%s" % (button_out_prefix,fn))
        button_out = Button(frame, text="%s%s" % (button_out_prefix,button_out_nofile), command=find_filename_out)
        button_out.pack(padx=3, pady=3)

        # handle run toggle
        check_run_var = IntVar(frame)
        check_run = Checkbutton(frame, text=HELP_TEXT_RUN, variable=check_run_var, onvalue=1, offvalue=0)
        check_run.pack()

        # handle AWS_ACCESS_KEY_ID
        entry_AWS_ACCESS_KEY_ID_default = "Enter AWS_ACCESS_KEY_ID (only necessary for 'reflow run')"
        entry_AWS_ACCESS_KEY_ID = Entry(frame, width=60)
        entry_AWS_ACCESS_KEY_ID.insert(END, entry_AWS_ACCESS_KEY_ID_default)
        entry_AWS_ACCESS_KEY_ID.pack()

        # handle AWS_SECRET_ACCESS_KEY
        entry_AWS_SECRET_ACCESS_KEY_default = "Enter AWS_SECRET_ACCESS_KEY (only necessary for 'reflow run')"
        entry_AWS_SECRET_ACCESS_KEY = Entry(frame, width=60)
        entry_AWS_SECRET_ACCESS_KEY.insert(END, entry_AWS_SECRET_ACCESS_KEY_default)
        entry_AWS_SECRET_ACCESS_KEY.pack()

        # handle AWS_REGION
        entry_AWS_REGION_default = "Enter AWS_REGION (only necessary for 'reflow run')"
        entry_AWS_REGION = Entry(frame, width=60)
        entry_AWS_REGION.insert(END, entry_AWS_REGION_default)
        entry_AWS_REGION.pack()

        # handle input RF file selection
        def find_filenames_rfs():
            rfs = askopenfilenames(title="Select Input RF File(s)", filetypes=(('RF file','*.rf'),))
            text_rfs.config(state='normal')
            text_rfs.delete('1.0', END)
            if len(rfs) == 0:
                text_rfs.insert('1.0', '<no RF file(s) selected>')
            else:
                text_rfs.insert('1.0', '\n'.join(rfs))
            text_rfs.config(state='disabled')
        button_rfs = Button(frame, text="Input RF file(s):", command=find_filenames_rfs)
        button_rfs.pack(padx=3, pady=3)
        text_rfs = ScrolledText(frame, height=10)
        text_rfs.insert('1.0', '<no RF file(s) selected>')
        text_rfs.config(state='disabled')
        text_rfs.pack()

        # handle save button
        def finish_applet():
            valid = True
            # check output batch RF
            try:
                if button_out['text'] == "%s%s" % (button_out_prefix,button_out_nofile):
                    gui_popup("ERROR: Output Batch RF file not selected", title="ERROR"); valid = False
            except:
                pass
            # check input sample RF(s)
            try:
                if text_rfs.get('1.0', END).strip() == '<no RF file(s) selected>':
                    gui_popup("ERROR: Input RF file(s) not selected", title="ERROR"); valid = False
            except:
                pass
            if valid:
                out_fn = button_out['text'].split(':')[-1].strip()
                sys.argv.append('-o'); sys.argv.append(out_fn)
                if isfile(out_fn):
                    remove(out_fn)
                if check_run_var.get() == 1:
                    sys.argv.append('--run')
                sys.argv += [rf.strip() for rf in text_rfs.get('1.0', END).strip().splitlines()]
                if entry_AWS_ACCESS_KEY_ID.get().strip() != entry_AWS_ACCESS_KEY_ID_default.strip():
                    environ['AWS_ACCESS_KEY_ID'] = entry_AWS_ACCESS_KEY_ID.get().strip()
                if entry_AWS_SECRET_ACCESS_KEY.get().strip() != entry_AWS_SECRET_ACCESS_KEY_default.strip():
                    environ['AWS_SECRET_ACCESS_KEY'] = entry_AWS_SECRET_ACCESS_KEY.get().strip()
                if entry_AWS_REGION.get().strip() != entry_AWS_REGION_default.strip():
                    environ['AWS_REGION'] = AWS_REGION.get().strip()
                try:
                    root.destroy()
                except:
                    pass
        button_save = Button(frame, text="Save Batch RF", command=finish_applet)
        button_save.pack(padx=3, pady=3)

        # add title and execute GUI
        root.title("Reflow Batch")
        root.mainloop()
    except:
        print("ERROR: Unable to import Tkinter", file=stderr); exit(1)
    if len(sys.argv) == 1:
        exit()

# run reflow as well
def run_reflow(batch_rf_fn):
    def execute(batch_rf_fn):
        # check if reflow is in PATH
        reflow_exe_path = None
        try:
            o,e = Popen(['reflow', 'run', '-h'], stdin=PIPE, stdout=PIPE, stderr=PIPE).communicate()
            if e.decode().startswith('usage: reflow run'):
                reflow_exe_path = 'reflow'
        except:
            pass

        # if reflow doesn't exist, set it up
        if reflow_exe_path is None:
            if GUI:
                text_stdout.insert(END, "Setting up reflow...\n")
            else:
                print("Setting up reflow...")
            if 'linux' in platform:
                url = REFLOW_EXEC_LINUX_URL
            elif 'darwin' in platform:
                url = REFLOW_EXEC_MAC_URL
            else:
                print("Your platform is not currently supported: %s" % platform)
            from urllib.request import urlopen
            reflow_dir = expanduser('~/.reflow')
            reflow_exe_path = '%s/reflow' % reflow_dir
            makedirs(reflow_dir, exist_ok=True)
            reflow_exe = urlopen(url).read()
            f = open(reflow_exe_path, 'wb'); f.write(reflow_exe); f.close()

        # run reflow
        command = [reflow_exe_path, 'run', batch_rf_fn]
        with Popen(command, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
            for line in p.stdout:
                if GUI:
                    text_stdout.insert(END, line)
                else:
                    print(line, end='')

    if GUI and False: # TODO REMOVE and False
        root = TK()
        root.geometry("600x400")
        frame = Frame(root)
        header = Label(frame, text="Reflow Console", font=('Arial',24))
        header.pack()
        text_stdout = ScrolledText(frame, height=10)
        text_stdout.pack()
        # TODO FINISH CONSOLE GUI
        execute(batch_rf_fn)
    else:
        execute(batch_rf_fn)

# main execution
if __name__ == "__main__":
    if len(sys.argv) == 1:
        GUI = True; run_gui()
    args = main()
    if args.run:
        run_reflow(args.output.name)
