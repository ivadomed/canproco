# Installation instructions

Installation

```
cd ~
mkdir anima/
cd anima/
wget -q https://github.com/Inria-Empenn/Anima-Public/releases/download/v4.2/Anima-macOS-4.2.zip # for MACOS
unzip Anima-macOS-4.2.zip
rm Anima-macOS-4.2.zip
git lfs install
git clone --depth 1 https://github.com/Inria-Visages/Anima-Scripts-Public.git
git clone --depth 1 https://github.com/Inria-Visages/Anima-Scripts-Data-Public.git
```

Configure directories

```
cd ~
mkdir .anima/
touch .anima/config.txt

echo "[anima-scripts]" >> .anima/config.txt
echo "anima = ${HOME}/anima/Anima-Binaries-4.2/" >> .anima/config.txt
echo "anima-scripts-public-root = ${HOME}/anima/Anima-Scripts-Public/" >> .anima/config.txt
echo "extra-data-root = ${HOME}/anima/Anima-Scripts-Data-Public/" >> .anima/config.txt
```
