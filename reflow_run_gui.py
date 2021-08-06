#! /usr/bin/env python3
'''
Simple GUI to help run 'reflow run'
'''
from os import chmod, environ, makedirs, stat
from os.path import expanduser
from stat import S_IEXEC
from subprocess import PIPE, Popen
from sys import platform, stderr

# helpful constants
REFLOW_EXEC_LINUX_URL = 'https://github.com/grailbio/reflow/releases/download/reflow1.6.0/reflow1.6.0.linux.amd64'
REFLOW_EXEC_MAC_URL = 'https://github.com/grailbio/reflow/releases/download/reflow1.6.0/reflow1.6.0.darwin.amd64'

# run GUI
def run_gui():
    try:
        # imports
        from tkinter import Button, Checkbutton, END, Entry, Frame, IntVar, Label, OptionMenu, StringVar, Tk
        from tkinter.filedialog import askopenfilename, asksaveasfilename
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
        header = Label(frame, text="Reflow Run GUI", font=('Arial',24))
        header.pack()

        # handle AWS_ACCESS_KEY_ID
        entry_AWS_ACCESS_KEY_ID_default = "Enter AWS_ACCESS_KEY_ID"
        entry_AWS_ACCESS_KEY_ID = Entry(frame, width=60)
        entry_AWS_ACCESS_KEY_ID.insert(END, entry_AWS_ACCESS_KEY_ID_default)
        entry_AWS_ACCESS_KEY_ID.pack()

        # handle AWS_SECRET_ACCESS_KEY
        entry_AWS_SECRET_ACCESS_KEY_default = "Enter AWS_SECRET_ACCESS_KEY"
        entry_AWS_SECRET_ACCESS_KEY = Entry(frame, width=60)
        entry_AWS_SECRET_ACCESS_KEY.insert(END, entry_AWS_SECRET_ACCESS_KEY_default)
        entry_AWS_SECRET_ACCESS_KEY.pack()

        # handle AWS_REGION
        entry_AWS_REGION_default = "Enter AWS_REGION"
        entry_AWS_REGION = Entry(frame, width=60)
        entry_AWS_REGION.insert(END, entry_AWS_REGION_default)
        entry_AWS_REGION.pack()

        # handle input RF file selection
        button_rf_prefix = "Input RF:\n"
        button_rf_nofile = "<none selected>"
        def find_filename_rf():
            fn = askopenfilename(title="Select Input RF File", filetypes=(('RF file','*.rf'),))
            if len(fn) == 0:
                button_rf.configure(text="%s%s" % (button_rf_prefix,button_rf_nofile))
            else:
                button_rf.configure(text="%s%s" % (button_rf_prefix,fn))
        button_rf = Button(frame, text="%s%s" % (button_rf_prefix,button_rf_nofile), command=find_filename_rf)
        button_rf.pack(padx=3, pady=3)

        # handle save button
        def finish_applet():
            valid = True
            # check input RF
            try:
                if button_rf['text'] == "%s%s" % (button_rf_prefix,button_rf_nofile):
                    gui_popup("ERROR: Input RF file not selected", title="ERROR"); valid = False
            except:
                pass
            # check AWS_ACCESS_KEY_ID
            try:
                if entry_AWS_ACCESS_KEY_ID.get().strip() == entry_AWS_ACCESS_KEY_ID_default.strip():
                    gui_popup("ERROR: AWS_ACCESS_KEY_ID not entered", title="ERROR"); valid = False
            except:
                pass
            # check AWS_SECRET_ACCESS_KEY
            try:
                if entry_AWS_SECRET_ACCESS_KEY.get().strip() == entry_AWS_SECRET_ACCESS_KEY_default.strip():
                    gui_popup("ERROR: AWS_SECRET_ACCESS_KEY not entered", title="ERROR"); valid = False
            except:
                pass
            # check AWS_REGION
            try:
                if entry_AWS_REGION.get().strip() == entry_AWS_REGION_default.strip():
                    gui_popup("ERROR: AWS_REGION not entered", title="ERROR"); valid = False
            except:
                pass
            if valid:
                # set environment variables
                text_log.configure(state='normal')
                text_log.insert(END, "Using AWS_ACCESS_KEY_ID: %s\n" % entry_AWS_ACCESS_KEY_ID.get().strip())
                environ['AWS_ACCESS_KEY_ID'] = entry_AWS_ACCESS_KEY_ID.get().strip()
                text_log.insert(END, "Using AWS_SECRET_ACCESS_KEY: %s\n" % entry_AWS_SECRET_ACCESS_KEY.get().strip())
                environ['AWS_SECRET_ACCESS_KEY'] = entry_AWS_SECRET_ACCESS_KEY.get().strip()
                text_log.insert(END, "Using AWS_REGION: %s\n" % entry_AWS_REGION.get().strip())
                environ['AWS_REGION'] = entry_AWS_REGION.get().strip()
                text_log.configure(state='disabled')
                root.update_idletasks()

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
                    text_log.configure(state='normal')
                    text_log.insert(END, "'reflow' executable not found in path\nDownloading reflow (may take a minute or so)...\n")
                    text_log.configure(state='disabled')
                    root.update_idletasks()
                    if 'linux' in platform:
                        url = REFLOW_EXEC_LINUX_URL
                    elif 'darwin' in platform:
                        url = REFLOW_EXEC_MAC_URL
                    else:
                        text_log.configure(state='normal')
                        text_log.insert(END, "Your platform is not currently supported: %s\n" % platform); return
                        text_log.configure(state='disabled')
                        root.update_idletasks()
                    from urllib.request import urlopen
                    reflow_dir = expanduser('~/.reflow')
                    reflow_exe_path = '%s/reflow' % reflow_dir
                    makedirs(reflow_dir, exist_ok=True)
                    reflow_exe = urlopen(url).read()
                    f = open(reflow_exe_path, 'wb'); f.write(reflow_exe); f.close()
                    chmod(reflow_exe_path, stat(reflow_exe_path).st_mode | S_IEXEC)
                    text_log.configure(state='normal')
                    text_log.insert(END, "Successfully downloaded reflow\nPerforming initial reflow configuration...\n")
                    root.update_idletasks()
                    return # TODO REMOVE WHEN DONE

                # end GUI
                try:
                    root.destroy()
                except:
                    pass
        button_run = Button(frame, text="Run", command=finish_applet)
        button_run.pack(padx=3, pady=3)

        # add console log
        header = Label(frame, text="\nLog:", font=('Arial',12))
        header.pack()
        text_log = ScrolledText(frame, height=10, state='disabled')
        text_log.pack()

        # add title and execute GUI
        root.title("Reflow Run GUI")
        root.mainloop()
    except Exception as e:
        print(e); exit()
        print("ERROR: Unable to import Tkinter", file=stderr); exit(1)

# main execution
if __name__ == "__main__":
    run_gui()
