#!/usr/bin/env bash

# used for travis only
if [[ -z "$TRAVIS_OS_NAME" ]]
then
    echo "[05-ssl-dev] only for use on travis ci!"
    exit 1
fi

if [[ $TRAVIS_OS_NAME == linux ]]
then
    sudo apt-get update
    sudo apt-get install -y libssl-dev
    echo "[05-ssl-dev] linux ssl headers installed"
    exit 0
fi

if [[ $TRAVIS_OS_NAME == osx ]]
then
    brew update
    brew install --force openssl
    echo "[05-ssl-dev] macos brew ssl headers installed"
    exit 0
fi

echo "[05-ssl-dev] unknown TRAVIS_OS_NAME '$TRAVIS_OS_NAME'"
exit 1