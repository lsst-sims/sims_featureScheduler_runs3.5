Don't forget to run the DDF batch with the old DDF positions.
UPDATE table SET new_column=old_column ;

sqlite3 early_ss_template_fpextra1.5_v3.4_10yrs.db "ALTER TABLE observations ADD scheduler_note VARCHAR;"  "UPDATE observations SET scheduler_note=note;" ".exit"