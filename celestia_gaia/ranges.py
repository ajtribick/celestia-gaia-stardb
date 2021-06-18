# gaia-stardb: Processing Gaia DR2 for celestia.Sci/Celestia
# Copyright (C) 2019–2021  Andrew Tribick
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""Range handling."""

from __future__ import annotations

from typing import Iterator


class Range:
    """Represents an inclusive integer range."""
    begin: int
    end: int

    def __init__(self, begin: int, end: int) -> None:
        """Creates a new Range."""
        assert begin <= end
        self.begin = begin
        self.end = end

    def __repr__(self) -> str:
        return f'Range({self.begin!r}, {self.end!r})'

    def __str__(self) -> str:
        return f'[{self.begin}, {self.end}]'

    def subtract_range(self, *args) -> list[Range]:
        """Returns the ranges that would result from removing a range from this one."""
        if len(args) == 1:
            if isinstance(args[0], Range):
                other = args[0]
            else:
                other = Range(args[0][0], args[0][1])
        elif len(args) == 2:
            other = Range(args[0], args[1])
        else:
            raise TypeError(f'subtract_range() takes 1 or 2 arguments ({len(args)} given)')

        if other.begin <= self.begin and other.end >= self.end:
            result = []
        elif other.end < self.begin or other.begin > self.end:
            result = [self]
        else:
            result = []
            if other.begin > self.begin:
                result.append(Range(self.begin, other.begin-1))
            if other.end < self.end:
                result.append(Range(other.end+1, self.end))
        return result

    def chunks(self, chunk_size) -> list[Range]:
        """Splits the range into chunks of size chunk_size, last chunk may be smaller."""
        result = []
        start = self.begin
        while start <= self.end:
            result.append(Range(start, min(start+chunk_size-1, self.end)))
            start += chunk_size
        return result


class MultiRange:
    """Represents a set of ranges."""
    ranges: list[Range]

    def __init__(self, begin: int, end: int) -> None:
        self.ranges = [Range(begin, end)]

    def __bool__(self) -> bool:
        return len(self.ranges) > 0

    def remove(self, *args) -> None:
        """Removes a range from the set of ranges."""
        if len(args) == 1:
            if isinstance(args[0], Range):
                other = args[0]
            else:
                other = Range(args[0][0], args[0][1])
        elif len(args) == 2:
            other = Range(args[0], args[1])
        else:
            raise TypeError(f'remove() takes 1 or 2 arguments ({len(args)} given)')

        new_ranges = []
        for subrange in self.ranges:
            new_ranges += subrange.subtract_range(other)

        self.ranges = new_ranges

    def chunk_ranges(self, chunk_size: int) -> Iterator[Range]:
        """Returns the constituent ranges with a maximum chunk size."""
        for subrange in self.ranges:
            for chunk in subrange.chunks(chunk_size):
                yield chunk
