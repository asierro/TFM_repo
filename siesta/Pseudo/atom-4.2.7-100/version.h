c     Version.h
c
c     This file MUST be updated after every self-consistent commit,
c     and the PL ("patch level") number increased by one, unless the
c     modification involves raising a minor or major version number,
c     in which case the PL should be reset to zero.
c
c     A self-consistent commit is a group of changes that fix a bug
c     or implement a new feature, in such a way that the program can
c     be compiled (no loose ends left). An update to the ChangeLog file
c     should be an integral part of a commit (the PL number should be
c     included for reference.)
c
c     After it is done, this file should be commited.
c
c     Please note that the use of a "data" statement forces the inclusion
c     of this file at the end of the declaration section, and before any
c     other data statements.
c
      character*60 version
      character*10 atom_id
      data version
     $     /"ATM Version 4.2.7 (revno 100)"/
      data atom_id /"ATM4.2.7"/

