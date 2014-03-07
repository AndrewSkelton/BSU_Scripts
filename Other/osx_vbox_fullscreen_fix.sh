
##OS X MAVERICKS VBOX FULLSCREEN SOLUTION
sudo defaults write /Applications/VirtualBox.app/Contents/Resources/vmstarter.app/Contents/Info.plist LSUIPresentationMode -int 4
sudo defaults write /Applications/VirtualBox.app/Contents/Resources/VirtualBoxVM.app/Contents/Info.plist LSUIPresentationMode -int 4
sudo chmod 644 /Applications/VirtualBox.app/Contents/Resources/VirtualBoxVM.app/Contents/Info.plist
sudo chmod 644 /Applications/VirtualBox.app/Contents/Resources/vmstarter.app/Contents/Info.plist