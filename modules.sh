#!/bin/bash
# Shared functions for scripts
# Copyright (C) 2022 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 12 Nov 2023; latest update: 12 Nov 2023

print_failure_message() {
    echo "Failure: assembly of isolate $1 could not be polished."
}
