#!/bin/bash

command_exists () {
    command -v $1 >/dev/null 2>&1
}

install_if_not_exists () {
  if ! command_exists $1 ; then
    echo "Installing $1..."
    sudo apt-get install "$1"
  fi
}

# Ruby
install_if_not_exists ruby
install_if_not_exists zenity

# Bundler
if ! gem list | grep bundler &> /dev/null ; then
  echo 'Installing bundler...'
  gem install bundler
fi

bundle install

echo 'Ready to run utils from "bin/" directory!'
