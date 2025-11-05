#!/usr/bin/env python3
"""
HDF5 Slice Viewer (2D/3D) with GUI
----------------------------------
- Pick an .h5/.hdf5 file via a file dialog
- Choose a dataset (label) from the tree
- If 2D -> show image (with colorbar)
- If 3D -> choose axis and move a slider to browse slices (colorbar updates)
- Axis labels in cm: x (cm), y (cm), z (cm)

Dependencies: tkinter, h5py, numpy, matplotlib
Run: python h5_3d_viewer.py
"""
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import numpy as np
import h5py
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

class H5SliceViewer(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("HDF5 2D/3D Viewer")
        self.geometry("1000x700")
        self.minsize(900, 600)

        # State
        self.h5_path = None
        self.h5_file = None
        self.datasets = []  # list of (path, shape, dtype)
        self.current_dataset_path = None
        self.current_data = None
        self.slice_axis = tk.StringVar(value="2")  # "0", "1", "2"
        self.slice_index = tk.IntVar(value=0)
        self.cbar = None  # persistent colorbar handle

        # UI
        self._build_ui()

        # Close hook
        self.protocol("WM_DELETE_WINDOW", self.on_close)

    # --------------- UI ----------------
    def _build_ui(self):
        # Top bar: open file + info
        top = ttk.Frame(self)
        top.pack(side=tk.TOP, fill=tk.X, padx=8, pady=8)

        open_btn = ttk.Button(top, text="Open HDF5 file…", command=self.open_file)
        open_btn.pack(side=tk.LEFT)

        self.file_lbl = ttk.Label(top, text="No file loaded", foreground="#444")
        self.file_lbl.pack(side=tk.LEFT, padx=10)

        # Main PanedWindow: left (datasets list) / right (plot + controls)
        pw = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        pw.pack(fill=tk.BOTH, expand=True)

        # Left: datasets tree
        left = ttk.Frame(pw)
        pw.add(left, weight=1)

        lbl = ttk.Label(left, text="Datasets")
        lbl.pack(anchor="w", padx=8, pady=(8, 0))

        self.tree = ttk.Treeview(left, columns=("shape","dtype"), show="tree headings", selectmode="browse")
        self.tree.heading("#0", text="Path")
        self.tree.heading("shape", text="Shape")
        self.tree.heading("dtype", text="DType")
        self.tree.column("#0", width=300, stretch=True)
        self.tree.column("shape", width=120, anchor="center")
        self.tree.column("dtype", width=100, anchor="center")
        self.tree.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        self.tree.bind("<<TreeviewSelect>>", self.on_tree_select)

        # Right: plot + controls
        right = ttk.Frame(pw)
        pw.add(right, weight=3)

        # Matplotlib Figure
        fig_frame = ttk.Frame(right)
        fig_frame.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        self.fig = Figure(figsize=(5,4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title("No dataset selected")

        self.canvas = FigureCanvasTkAgg(self.fig, master=fig_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Controls for 3D slicing
        ctrl = ttk.LabelFrame(right, text="Controls")
        ctrl.pack(fill=tk.X, padx=8, pady=(0,8))

        # Axis selection
        ttk.Label(ctrl, text="Slice axis:").grid(row=0, column=0, padx=6, pady=6, sticky="w")
        axis_combo = ttk.Combobox(ctrl, state="readonly", textvariable=self.slice_axis, values=["0","1","2"], width=5)
        axis_combo.grid(row=0, column=1, padx=6, pady=6, sticky="w")
        axis_combo.bind("<<ComboboxSelected>>", lambda e: self.update_slice_from_controls())

        # Slider
        ttk.Label(ctrl, text="Slice index:").grid(row=0, column=2, padx=6, pady=6, sticky="w")
        self.slice_slider = ttk.Scale(ctrl, from_=0, to=0, orient=tk.HORIZONTAL, command=self.on_slider_move, length=300)
        self.slice_slider.grid(row=0, column=3, padx=6, pady=6, sticky="we")

        # Index readout
        self.slice_idx_lbl = ttk.Label(ctrl, text="0")
        self.slice_idx_lbl.grid(row=0, column=4, padx=6, pady=6, sticky="w")

        ctrl.columnconfigure(3, weight=1)

        # Help footer
        footer = ttk.Frame(self)
        footer.pack(side=tk.BOTTOM, fill=tk.X, padx=8, pady=(0,8))
        ttk.Label(footer, text="Tips: 2D datasets are displayed as images. For 3D datasets, pick an axis and move the slider to browse slices.").pack(side=tk.LEFT)

    # ------------- HDF5 handling --------------
    def open_file(self):
        path = filedialog.askopenfilename(
            title="Select HDF5 file",
            filetypes=[("HDF5 files", "*.h5 *.hdf5"), ("All files", "*.*")]
        )
        if not path:
            return
        if self.h5_file is not None:
            try:
                self.h5_file.close()
            except Exception:
                pass
            self.h5_file = None

        try:
            self.h5_file = h5py.File(path, "r")
        except Exception as e:
            messagebox.showerror("Error opening file", f"Failed to open file:\n{e}")
            return

        self.h5_path = path
        self.file_lbl.config(text=os.path.basename(path))
        self.populate_datasets()
        self.current_dataset_path = None
        self.current_data = None
        self.ax.clear()
        self.ax.set_title("Select a dataset from the left")
        self._clear_colorbar()
        self.canvas.draw()

    def populate_datasets(self):
        # Clear tree
        for i in self.tree.get_children():
            self.tree.delete(i)
        self.datasets.clear()
        if self.h5_file is None:
            return

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                shape = obj.shape
                dtype = str(obj.dtype)
                self.datasets.append((name, shape, dtype))
        self.h5_file.visititems(visit)

        # Build tree groups by path
        groups = {}
        for path, shape, dtype in self.datasets:
            parts = path.strip("/").split("/")
            parent_id = ""
            for i, part in enumerate(parts):
                subpath = "/".join(parts[:i+1])
                if i < len(parts) - 1:
                    if subpath not in groups:
                        parent_tree_id = "" if i == 0 else groups["/".join(parts[:i])]
                        groups[subpath] = self.tree.insert(parent_tree_id, "end", text=part, values=("", ""))
                else:
                    parent_tree_id = "" if i == 0 else groups["/".join(parts[:i])]
                    self.tree.insert(parent_tree_id, "end", text=part, values=(str(shape), dtype), iid=path)

        if not self.datasets:
            messagebox.showinfo("No datasets", "No datasets found in this file.")

    def on_tree_select(self, event):
        sel = self.tree.selection()
        if not sel:
            return
        path = sel[0]  # iid is dataset path
        vals = self.tree.item(path, "values")
        if not vals:
            return  # clicked a group node
        self.load_dataset(path)

    def load_dataset(self, dset_path):
        if self.h5_file is None:
            return
        try:
            dset = self.h5_file[dset_path]
            data = dset[()]
        except Exception as e:
            messagebox.showerror("Error reading dataset", f"Failed to read dataset:\n{e}")
            return

        if not isinstance(data, np.ndarray):
            messagebox.showwarning("Unsupported", "Selected dataset is not a numeric array.")
            return

        self.current_dataset_path = dset_path
        self.current_data = np.asarray(data)

        # Remove leading singleton dimensions (except keep 2D/3D shape semantics)
        while self.current_data.ndim > 2 and self.current_data.shape[0] == 1:
            self.current_data = self.current_data[0]

        self.update_plot_initial()

    # -------------- Plotting -----------------
    def _clear_colorbar(self):
        if self.cbar is not None:
            try:
                self.cbar.remove()
            except Exception:
                pass
            self.cbar = None

    def _update_colorbar(self, mappable):
        # Remove old colorbar if exists, then create a new one
        self._clear_colorbar()
        try:
            self.cbar = self.fig.colorbar(mappable, ax=self.ax, fraction=0.046, pad=0.04)
        except Exception:
            self.cbar = None

    def update_plot_initial(self):
        if self.current_data is None:
            return
        arr = self.current_data
        self.ax.clear()
        self._clear_colorbar()

        if arr.ndim == 1:
            self.ax.plot(arr)
            self.ax.set_title(f"{self.current_dataset_path} (1D)")
            self._set_controls_enabled(False)
            self.ax.set_xlabel("index")
            self.ax.set_ylabel(self.current_dataset_path.split("/")[-1])
        elif arr.ndim == 2:
            im = self.ax.imshow(arr, origin="lower", aspect="auto", cmap="viridis")
            self._update_colorbar(im)
            self.ax.set_title(f"{self.current_dataset_path} (2D)")
            self._set_controls_enabled(False)
            self.ax.set_xlabel("x (cm)")
            self.ax.set_ylabel("y (cm)")
        elif arr.ndim == 3:
            # Default: slice along axis 2 (k)
            self.slice_axis.set("2")
            kmax = max(0, arr.shape[2]-1)
            self.slice_slider.configure(from_=0, to=kmax)
            self.slice_index.set(int(kmax//2))
            self.slice_idx_lbl.config(text=str(self.slice_index.get()))
            self._set_controls_enabled(True)
            self.draw_slice()
            return  # draw_slice handles colorbar and labels
        else:
            self.ax.text(0.5, 0.5, f"ndim={arr.ndim} not supported", ha="center", va="center", transform=self.ax.transAxes)
            self._set_controls_enabled(False)

        self.canvas.draw()

    def _set_controls_enabled(self, enabled: bool):
        # Enable/disable the slider. Axis combo stays usable only for 3D;
        # simplest approach: disable slider when not 3D.
        if enabled:
            self.slice_slider.state(["!disabled"])
        else:
            self.slice_slider.state(["disabled"])

    def on_slider_move(self, val):
        try:
            idx = int(float(val))
        except Exception:
            idx = 0
        self.slice_index.set(idx)
        self.slice_idx_lbl.config(text=str(idx))
        self.draw_slice()

    def update_slice_from_controls(self):
        if self.current_data is None or self.current_data.ndim != 3:
            return
        axis = int(self.slice_axis.get())
        n = self.current_data.shape[axis]
        self.slice_slider.configure(from_=0, to=max(0, n-1))
        if self.slice_index.get() >= n:
            self.slice_index.set(n//2)
        self.slice_idx_lbl.config(text=str(self.slice_index.get()))
        self.draw_slice()

    def draw_slice(self):
        if self.current_data is None or self.current_data.ndim != 3:
            return
        arr = self.current_data
        axis = int(self.slice_axis.get())
        idx = self.slice_index.get()

        # Bounds check
        n = arr.shape[axis]
        if n == 0:
            return
        if idx >= n:
            idx = n - 1
            self.slice_index.set(idx)
            self.slice_idx_lbl.config(text=str(idx))

        # Extract slice
        if axis == 0:
            sl = arr[idx, :, :]
        elif axis == 1:
            sl = arr[:, idx, :]
        else:
            sl = arr[:, :, idx]

        self.ax.clear()
        im = self.ax.imshow(sl, origin="lower", aspect="auto", cmap="viridis")
        self._update_colorbar(im)
        self.ax.set_title(f"{self.current_dataset_path} (3D) axis={axis} idx={idx}")

        # Dynamic axis labels based on slicing axis
        if axis == 0:  # slicing X
            self.ax.set_xlabel("y (cm)")
            self.ax.set_ylabel("z (cm)")
        elif axis == 1:  # slicing Y
            self.ax.set_xlabel("x (cm)")
            self.ax.set_ylabel("z (cm)")
        else:  # slicing Z
            self.ax.set_xlabel("x (cm)")
            self.ax.set_ylabel("y (cm)")

        self.canvas.draw()

    # -------------- Cleanup ------------------
    def on_close(self):
        try:
            if self.h5_file is not None:
                self.h5_file.close()
        except Exception:
            pass
        self.destroy()

def main():
    app = H5SliceViewer()
    app.mainloop()

if __name__ == "__main__":
    main()
