#! /usr/bin/env python3
'''
Create a Reflow run file that batch-executes many other Reflow run files.
'''
from os import remove
from os.path import abspath, isfile
from sys import stderr, stdout
import argparse
import sys

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
        root.geometry("600x600")

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
        text_rfs = ScrolledText(frame, height=5)
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
    print(batch_rf_fn)
    exit(1) # TODO

# main execution
if __name__ == "__main__":
    if len(sys.argv) == 1:
        run_gui()
    args = main()
    if args.run:
        run_reflow(args.output.name)
