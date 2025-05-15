# INSTALLING


## Step 1: latticelab, a required dependency
```bash
git clone https://github.com/Spuriosity1/latticelab.git && cd latticelab
meson setup build 
ninja -C build install
```
If you can't touch `/usr/local`, instead do
```bash
meson setup build -Dprefix="/your/install/prefix"
```

## Step 2: install this

```bash
git clone https://github.com/Spuriosity1/dmnd_dilute.git && cd dmnd_dilute
meson setup build
# if this fails to find latticelab, rune 
export PKG_CONFIG_PATH="/your/install/prefix/lib/pkgconfig"
ninja -C build
```
