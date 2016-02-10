;#####################################################################
; Copyright 2004-2006.
;#####################################################################
; physbam.el
;#####################################################################

;#####################################################################
; Bind Projects to Viewers
;#####################################################################
(setq physbam-projects-to-viewers
      '(("solids_3d" . "opengl_3d")
        ("solids_2d" . "opengl_2d")
        ("levelset_quadtree" . "opengl_2d")
        ("fluids_2d" . "opengl_2d")
        ("fluids_quadtree" . "opengl_2d")
        ("fluids_3d" . "opengl_3d")
        ("fluids_octree" . "opengl_3d")
        ("smoke_and_fire_2d" . "opengl_2d")
        ("smoke_and_fire_quadtree" . "opengl_2d")
        ("smoke_and_fire_3d" . "opengl_3d")
        ("smoke_and_fire_octree" . "opengl_3d")
        ("water_free_surface_2d" . "opengl_2d")
        ("water_free_surface_quadtree" . "opengl_2d")
        ("water_free_surface_3d" . "opengl_3d")
        ("face" . "opengl_3d")
        ("articulated_rigid_bodies" . "opengl_3d")
        ("water_free_surface_octree" . "opengl_3d")))

;#####################################################################
; Compiler Table
;#####################################################################
(setq physbam-compilers
      '(("gcc-3.4" . ("/usr/local/compilers/gcc-3.4/bin/gcc" 
                      "/usr/local/compilers/gcc-3.4/bin/g++" 
                      "/usr/local/compilers/icecream/gcc-3.4.tar.bz2"))
        ("gcc-3.4-64" .  ("/usr/local/compilers/gcc-3.4-64/bin/gcc" 
                          "/usr/local/compilers/gcc-3.4-64/bin/g++" 
                          "/usr/local/compilers/icecream/gcc-3.4-64.tar.bz2,i686:/usr/local/compilers/icecream/gcc-3.4-64-cross.tar.bz2"))
        ("gcc-3.3.2" .  ("/usr/local/compilers/gcc-3.3.2/bin/gcc" 
                         "/usr/local/compilers/gcc-3.3.2/bin/g++" 
                         "/usr/local/compilers/icecream/gcc-3.3.2.tar.bz2"))
        ("gcc-4.0.1" .  ("/usr/local/compilers/gcc-4.0.1-i686-i686/bin/gcc" 
                         "/usr/local/compilers/gcc-4.0.1-i686-i686/bin/g++" 
                         "/usr/local/compilers/icecream/gcc-4.0.1-i686-i686.tar.bz2"))
        ("gcc-4.1.1" .  ("/usr/local/compilers/gcc-4.1.1-i686-i686/bin/gcc" 
                         "/usr/local/compilers/gcc-4.1.1-i686-i686/bin/g++" 
                         "/usr/local/compilers/icecream/gcc-4.1.1-i686-i686.tar.bz2"))
        ("gcc-4.0.1-64" .  ("/usr/local/compilers/gcc-4.0.1-x86_64-x86_64/bin/gcc" 
                            "/usr/local/compilers/gcc-4.0.1-x86_64-x86_64/bin/g++" 
                            "/usr/local/compilers/icecream/gcc-4.0.1-x86_64-x86_64.tar.bz2,i686:/usr/local/compilers/icecream/gcc-4.0.1-i686-x86_64.tar.bz2"))
        ("gcc-4.1.1-64" .  ("/usr/local/compilers/gcc-4.1.1-x86_64-x86_64/bin/gcc" 
                            "/usr/local/compilers/gcc-4.1.1-x86_64-x86_64/bin/g++" 
                            "/usr/local/compilers/icecream/gcc-4.1.1-x86_64-x86_64.tar.bz2,i686:/usr/local/compilers/icecream/gcc-4.1.1-i686-x86_64.tar.bz2"))
        ("icc" .  ("icc" "icc" "none"))))

;#####################################################################
; PhysBAM Style
;#####################################################################

(defconst physbam-c-style
  '((c-basic-offset . 4)
        (c-offsets-alist . ((string . -1000)
                        (c . c-lineup-C-comments)
                        (defun-open . 0)
                        (defun-close . 0)
                        (defun-block-intro . +)
                        (class-open . 0)
                        (class-close . 0)
                        (inline-open . 0)
                        (inline-close . 0)
                        (func-decl-cont . +)
                        (knr-argdecl-intro . 5)
                        (knr-argdecl . 0)
                        (topmost-intro . 0)
                        (topmost-intro-cont . 0)
                        (member-init-intro . +)
                        (member-init-cont . -1)
                        (inher-intro . +)
                        (inher-cont . c-lineup-multi-inher)
                        (block-open . 0)
                        (block-close . 0)
                        (brace-list-open . 0)
                        (brace-list-close . 0)
                        (brace-list-intro . +)
                        (brace-list-entry . 0)
                        (statement . 0)
                        (statement-cont . +)
                        (statement-block-intro . +)
                        (statement-case-intro . +)
                        (statement-case-open . +)
                        (substatement . +)
                        (substatement-open . 0)
                        (case-label . +)
                        (access-label . -)
                        (label . *)
                        (do-while-closure . 0)
                        (else-clause . 0)
                        (comment-intro . c-lineup-comment)
                        (arglist-intro . +)
                        (arglist-cont . 0)
                        (arglist-cont-nonempty . +)
                        (arglist-close . 0)
                        (stream-op . c-lineup-streamop)
                        (inclass . +)
                        (cpp-macro . -1000)
                        (friend . 0)
                        (extern-lang-open . 0)
                        (extern-lang-close . 0)
                        (inextern-lang . +)
			(namespace-open . 0)
			(namespace-close . 0)
			(innamespace . 0)
                        (template-args-cont . +)))    
    (c-hanging-braces-alist . ((class-open before after)
                               (class-close before)
                               (defun-open before after)
                               (defun-close before after)
                               (inline-open before)
                               (inline-close after)
                               (brace-list-open)
                               (brace-list-close)
                               (brace-list-intro)
                               (brace-list-entry)
                               (block-open after)
                               (block-close)
                               (substatement-open after)
                               (substatement-case-open)
			       (namespace-open after)
			       (namespace-close before after)
                               (extern-lang-open after)
                               (extern-lang-close after)))
    (c-hanging-colons-alist . ((case-label)
                               (label)
                               (access-label after)
                               (member-init-intro before)
                               (inher-intro)))
    (c-cleanup-list         . ((defun-close-semi)
                               (scope-operator)))
    (c-hanging-semi&comma-criteria . (c-semi&comma-no-newlines-before-nonblanks))
    (c-echo-syntactic-information-p . t))
  "PhysBAM C++ Style")

(c-add-style "physbam" physbam-c-style)

;#####################################################################
; PhysBAM Syntax Highlighting
;#####################################################################

(setq font-lock-keyword-case-fold-search nil) ; Need to be case sensitive
(font-lock-add-keywords 'c++-mode  '(("[^[:lower:]]\\([[:upper:]][[:upper:][:digit:]_]*\\)[<> ,]" 1 font-lock-type-face t))) ; match class names

(make-face 'font-lock-operators-face)
(font-lock-add-keywords 'c++-mode '(("[;<>:=!]\\|->" . 'font-lock-operators-face)))

(make-face 'font-lock-preprocessor-face)
(font-lock-add-keywords 'c++-mode '(("\\#[a-zA-Z0-9]*[ ]" . 'font-lock-preprocessor-face)))


;#####################################################################
; PhysBAM Helper Routines
;#####################################################################

(defun physbam-reduce (f x)
  (if (eq (cdr x) nil)
      (car x)
    (funcall f (car x) (physbam-reduce f (cdr x))))) 

(defun physbam-filter (f list)
  (let (filtered)
    (dolist (x list filtered)
      (if (funcall f x) (setq filtered (cons x filtered)) nil))))

;#####################################################################
; PhysBAM Navigation Commands
;#####################################################################

(defun physbam-header-flip ()
  "Find the .h file for this .C file (or vice versa)."
  (interactive)
  (let ((dotc (string-match "[.]\\(cpp\\|c\\)$" (buffer-file-name)))
        (doth (string-match "[.]h$" (buffer-file-name))))
    (if dotc
        (find-file (concat (substring (buffer-file-name) 0 dotc) ".h"))
      (if doth
          (find-file (concat (substring (buffer-file-name) 0 doth) ".cpp"))
        (message "Not a cpp or h file!!")))))

(defun physbam-dimension-flip ()
  "Find the .h file for this .C file (or vice versa)."
  (interactive)
  (let ((buf (buffer-name)))
    (if (string-match "\\(.+\\)\\(1D\\|2D\\|3D\\)\\(.+\\)"  buf)
      (let ((current (match-string 2 buf))
            (left (match-string 1 buf))
            (right (match-string 3 buf)))
        (find-file (concat left
                (cond ((string= "1D" current) "2D") ((string= "2D" current) "3D") ((string= "3D" current) "1D"))
                right))))))

(defun physbam-update-tags ()
  "Run update-tags.sh script to rebuild"
  (interactive) 
  (call-process "~/bin/update_tags.sh" nil nil)
  (kill-buffer "TAGS"))


(defun physbam-setup-filename-search ()
  (interactive)
  (if (not (eq (get-buffer "*physbam-filenames*") nil))
      (kill-buffer "*physbam-filenames*"))
  (call-process "find" nil "*physbam-filenames*" nil
                (getenv "PHYSBAM") "-name" "*.h" "-or" "-name" "*.cpp")
  (message "Done"))

(defun physbam-open-parent ()
  "Open header of parent class"
  (interactive)
  (let ((doth (string-match "[.]h$" (buffer-file-name))))
    (if doth
        (let ((save_point (point)))
          (goto-char 0)
          (if (re-search-forward ":\\(?:public\\|protected\\|private\\) \\(\\(?:\\w\\|_\\)*\\)" nil t)
              (let ((filename (concat (match-string 1) ".h")))
                (goto-char save_point)
                (if (file-exists-p filename)
                    (find-file filename)
                  (let ((physbam_filename (find-physbam-file filename)))
                    (if (file-exists-p physbam_filename)
                        (find-file physbam_filename)
                      (message "Can't find header of parent")))))
            (goto-char save_point)
            (message "No parent class")))
      (message "Not a header file"))))

(defun find-physbam-file (filename)
  (call-process "find" nil "*physbam-filename*" nil  (getenv "PUBLIC") "-name" filename)
  (set-buffer "*physbam-filename*")
  (goto-char 0)
  (let ((found_filename (buffer-substring (point) (line-end-position))))
    (kill-buffer "*physbam-filename*")
    (if (file-exists-p found_filename) found_filename nil)))

(defun physbam-grep-cpp-and-headers (directory_prefix querystr)
  (let ((path (concat (getenv "PHYSBAM") "/" directory_prefix "/")))
  (grep (concat "(find " path " -name '*.h' -o -name '*.cpp' | xargs grep -n -e \"" querystr "\") # "))))

(defun physbam-grep-public (querystr)
  "grep in Public_Library"
  (interactive "sGrep $PHYSBAM/Public_Library for:")
  (physbam-grep-cpp-and-headers "Public_Library" querystr))

(defun physbam-grep-projects (querystr)
  "grep in Projects"
  (interactive "sGrep $PHYSBAM/Projects for:")
  (physbam-grep-cpp-and-headers "*Projects" querystr))

(defun physbam-grep-tools (querystr)
  "grep in Tools"
  (interactive "sGrep $PHYSBAM/Tools for:")
  (physbam-grep-cpp-and-headers "Tools" querystr))

(defun physbam-fix-includes ()
  "Fix includes in current file"
  (interactive)
  (save-buffer)
  (shell-command (format "%s/Scripts/misc/fix_headers_file.sh %s" (getenv "PHYSBAM") (buffer-file-name)))
  (revert-buffer-no-prompt))

;#####################################################################
; PhysBAM Formatting Helpers
;#####################################################################

(defun physbam-fix-function-comment ()
  "Fix the physbam style function comments"
  (interactive)
  (setq function-name (physbam-get-function-name))
  (setq function-type (physbam-get-function-type function-name))
  (let ((line-count 0))
    (beginning-of-line)
    (if (looking-at "[a-zA-Z0-9_]") 
        (progn (forward-line -1)
               (setq line-count (+ line-count 1))))
    (while (looking-at "//") 
      (if (looking-at "//#")
          (while (looking-at "//")
            (delete-region (point) (line-end-position))
            (delete-char 1)
            (forward-line -1))
        (progn
          (setq line-count (+ line-count 1))
          (forward-line -1))))
    (forward-line +1)
    (physbam-insert-function-comment function-name function-type)
    (forward-line line-count)))

(defun physbam-insert-function-comment (name type)
  "Inserts a physbam style function comment at point"
  (interactive)
  (insert "//#####################################################################\n")
  (insert (concat "// " type))
  (if (string= type "Function")
      (insert (concat " " name)))
  (insert "\n")
  (insert "//#####################################################################\n"))

(defun physbam-get-function-name ()
  "Gets the name of the function on the current line"
  (beginning-of-line)
  (let ((curr-line (buffer-substring (point) (line-end-position))))
    (string-match "^\\(~?[a-zA-Z0-9_]+\\)\(.+$" curr-line)
    (match-string 1 curr-line)))

(defun physbam-get-function-type (name)
  "Determines the type of function on the current line"
  (forward-line -1)
  (beginning-of-line)
  (let ((curr-line (buffer-substring (point) (line-end-position))))
    (cond ((string-match "~" name) "Destructor")
          ((string-match name curr-line) "Constructor")
          (t "Function"))))

(defun physbam-get-lastname (fullname)
  (car (last (split-string fullname "[ ]+"))))

(defun physbam-current-year ()
  "Get year"
  (string-to-number (substring (current-time-string) -4)))

(defun physbam-fix-copyright-user (fullname)
  "Fix the list of names in copyright i.e. insert mine"
  (interactive)
  (goto-line 2)
  (beginning-of-line)
  (let ((first-line (buffer-substring (point) (line-end-position))))
    (let ((items (split-string first-line "[ \t]*,[ \t]*")))
      (let ((names '()) (years '()) (year-regexp "\\([0-9]\\{4\\}\\)"))
        (dolist (x items)
          (if (string-match (format "^\\(// Copyright[ ]\\)?+%s$" year-regexp) x)
              (setq years (cons (match-string 2 x) years))
            (if (string-match (format "^\\(// Copyright[ ]\\)?+%s-%s$" year-regexp year-regexp) x)
                (setq years (cons (match-string 2 x) (cons (match-string 3 x) years)))
              (setq names (cons (replace-regexp-in-string  "\\." "" x) names)))))
        (unless (member fullname names) (setq names (cons fullname names)))
        (let ((sorted-names (sort names '(lambda (x y) (string< (physbam-get-lastname x) (physbam-get-lastname y))))))
          (delete-region (point) (line-end-position))
          (insert (format "// Copyright %s, %s." 
                          (let ((min-year (physbam-reduce (lambda (x y) (if (< x y) x y)) (mapcar 'string-to-number years)))
                                (max-year (physbam-current-year)))
                            (if (= min-year max-year) (format "%d" min-year) (format "%d-%d" min-year max-year)))
                          (physbam-reduce (lambda (x y) (concat x ", " y)) sorted-names))))))))

(defun physbam-fix-copyright ()
  (interactive)
  (physbam-fix-copyright-user user-full-name))

(defun physbam-fix-copyright-user-list ()
  (interactive)
  (dolist (fullname user-list)
    (physbam-fix-copyright-user fullname)))

(defun physbam-insert-copyright ()
  (interactive)
  (goto-char 0) (goto-line 7)
  (insert "//#####################################################################\n")
  (insert (concat "// Copyright " (number-to-string (physbam-current-year))  ", " user-full-name ".\n"))
  (insert "// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.\n")
  (insert "//#####################################################################\n"))

(defun physbam-insert-header ()
  (interactive)
  (physbam-insert-copyright)
  (let ((classname (substring (buffer-name) 0 -2)))
    (insert (format "#ifndef __%s__\n" classname))
    (insert (format "#define __%s__\n\n" classname))
    (insert "namespace PhysBAM{\n\n")
    (insert "template<class T>\n")
    (insert (format "class %s\n" classname))
    (insert "{\n")
    (insert "public:\n\n")
    (insert "//#####################################################################\n")
    (insert "};\n")
    (insert "}\n")
    (insert "#endif\n")))

(defun physbam-add-include ()
  "Add include from current location"
  (interactive)
  (search-forward-regexp "[^A-Za-z0-9_]")
  (backward-char)
  (setq a (point))
  (search-backward-regexp "[^A-Za-z0-9_]")
  (forward-char)
  (setq b (point))
  (setq class-name (buffer-substring a b))
  (goto-line 0)
  (beginning-of-line)
  (search-forward "#include")
  (beginning-of-line)
  (insert (concat "#include <temp/" class-name ".h>\n"))
  (physbam-fix-includes))

(defun physbam-make-cpp-from-header ()
  (interactive)
  (goto-char 0)
  (goto-line 7)
  (copy-region-as-kill (point-min) (point))
  (physbam-header-flip)
  (yank)
  (let ((classname (substring (buffer-name) 0 -4)))
    (insert (format "#include \"%s.h\"\n" classname))
    (insert "using namespace PhysBAM;\n")
    (insert (format "template class %s<float>;\n" classname))
    (insert (format "template class %s<double>;\n" classname))))

(defun physbam-shift-indent-left ()
  (interactive) 
  (indent-rigidly (region-beginning) (region-end) -4))

(defun physbam-shift-indent-right ()
  (interactive) 
  (indent-rigidly (region-beginning) (region-end) 4))

;#####################################################################
; PhysBAM Build / Run Commands
;#####################################################################
; returns the last eleement of the path i.e. the project
(defun physbam-project-name (directory)
  (car (last (split-string directory "/" t))))

; Returns executable name in form <proj>_<release/debug>_<PLATFORM>  
(defun physbam-executable-name (base-executable use-release-always)
  (concat base-executable
          (if (or (string= (getenv "PLATFORM") nil) (string= (getenv "PLATFORM") "pentium4")) "" (concat "_" (getenv "PLATFORM")))
          (if (or use-release-always (string= physbam-project-type "release")) "" (concat "_" physbam-project-type))))

; runs a simulation
(setq-default physbam-run-parameters-history nil)
(defun physbam-run (run-params)
  (interactive (list (let ((execname (physbam-executable-name (physbam-project-name physbam-project-directory) nil))
                           (last-command (if physbam-run-parameters-history (car physbam-run-parameters-history) nil)))
                       (read-string (format "Run %s with%s as sim %s: "
                                            execname (if last-command (format " (%s)" last-command) "") current-prefix-arg)
                                    (if last-command last-command "")
                                    'physbam-run-parameters-history))))
  (let ((sim-number current-prefix-arg))
    (let* ((execname (physbam-executable-name (physbam-project-name physbam-project-directory) nil))
           (command (if physbam-tee-output
                        (format "cd %s; ./%s %s | tee output-%s.txt &"  physbam-project-directory execname run-params sim-number)
                      (format "cd %s; ./%s %s > output-%s.txt &"  physbam-project-directory execname run-params sim-number))))
      (shell-command command (concat "simulation-" (physbam-project-name physbam-project-directory) (if sim-number (format "-%s" sim-number) ""))))))

(defun physbam-run-debug ()
  (interactive)
  (let ((execname (physbam-executable-name (physbam-project-name physbam-project-directory) nil))
        (cmd (format "gdb --annotate=3 -d %s -cd=%s %s/%s" physbam-project-directory physbam-project-directory physbam-project-directory (physbam-executable-name (physbam-project-name physbam-project-directory) nil))))
    (message cmd)
    (gdb cmd)))

; runs the viewer
(defun physbam-run-viewer ()
  (interactive)
  (if (eq nil (assoc (physbam-project-name  physbam-project-directory) physbam-projects-to-viewers))
      (message (format "No viewer associated with project %s"
                       physbam-project-directory))
    (let ((viewer-name (cdr (assoc (physbam-project-name  physbam-project-directory) physbam-projects-to-viewers))))
      (shell-command (format "cd %s; %s . &" physbam-output-directory
                             (format "%s/Projects/%s/%s" (getenv "PHYSBAM") viewer-name (physbam-executable-name viewer-name t)))))))

(defun physbam-set-project-type (type)
  (setq physbam-project-type type)
  (physbam-setup-compile-command t))

(defun physbam-set-compile-mode (val)
  (setq physbam-compile-mode val)
  (physbam-setup-compile-command nil))

(defun physbam-set-project-directory (project-directory)
  (setq physbam-project-directory project-directory)
  (setq physbam-project-directory (concat physbam-project-directory (if (string= (substring physbam-project-directory -1 nil) "/") "" "/")))
  (physbam-read-project-settings)
  (physbam-setup-compile-command nil))

(defun physbam-choose-project-directory ()
  (interactive)
  (physbam-set-project-directory (read-file-name "New Project Directory: " physbam-project-directory () nil)))
  

(defun physbam-choose-output-directory ()
  (interactive)
  (setq physbam-output-directory (read-file-name "New output directory: " physbam-output-directory () nil))
  (physbam-setup-compile-command t))

(defun physbam-update-output-directory-menu ()
  (let* ((output-list (physbam-filter (lambda (x) 
                                        (and 
                                         (file-directory-p x) 
                                         (not (string= "." (substring x -1 nil))) 
                                         (not (string= ".." (substring x -2 nil)))
                                         (not (string= "CVS" (substring x -3 nil)))))
                                      (directory-files physbam-project-directory t)))
         (keymap-form (mapcar (lambda (x) `(,(file-name-nondirectory x) ,(file-name-nondirectory x) (nil) . (lambda () (interactive) (setq physbam-output-directory ,(concat x "/output")) (physbam-setup-compile-command t)))) output-list)))
    (define-key global-map [menu-bar physbam output-directory-list] 
      (cons '"Output Directories" (cons 'keymap (cons '"Select Output Directory" keymap-form))))))

(defun physbam-set-compile-count (count)
  (interactive)
  (setq physbam-compile-count count)
  (physbam-setup-compile-command t))

(defun physbam-set-compiler (compiler)
  (interactive)
  (setq physbam-compiler-id compiler)
  (let ((compiler-info (assoc physbam-compiler-id physbam-compilers)))
    (setq physbam-compiler (car (cdr (cdr compiler-info))))
    (setenv "ICECC_CXX" (car (cdr (cdr compiler-info))))
    (setenv "ICECC_CC" (car (cdr compiler-info)))
    (setenv "ICECC_VERSION" (car (cdr (cdr (cdr compiler-info))))))
  (physbam-setup-compile-command t))

(defun physbam-compile ()
  (interactive)
  (let ((old-default-directory default-directory))
    (setq default-directory physbam-project-directory)
    (call-interactively 'compile)
    (setq default-directory old-default-directory)))

(defun physbam-compile-current-file ()
  (interactive)
  (let ((compile-command (concat physbam-compiler " -I$PHYSBAM/Public_Library -DCOMPILE_WITHOUT_DYADIC_SUPPORT -DCOMPILE_WITHOUT_RLE_SUPPORT -c " (buffer-file-name))))
    (save-buffer)
    (call-interactively 'compile)
    (physbam-setup-compile-command nil)))

(defun physbam-compile-current-file-msvc ()
  (interactive)
  (let ((compile-command (concat "clwrap /nologo /MD /TP /O2 /GR /W3 /Wp64 /wd4996 /wd4355 /wd4150 /WL /EHsc /MD /DWIN32 /I$PHYSBAM/External_Libraries/boost /I$PHYSBAM/Public_Library /c " (file-name-nondirectory (buffer-file-name)) " /Fotest.obj " )))
    (save-buffer)
    (call-interactively 'compile)
    (physbam-setup-compile-command nil)))

(defun physbam-setup-compile-command (write_settings)
  (if physbam-use-scons
      (setq compile-command (format "nice scons --warn=no-duplicate-environment --warn=no-deprecated -Q --implicit-cache -k -u TYPE=%s %s"
                                physbam-project-type
                                (cond ((string= physbam-compile-mode "distcc") (format "CXX=\"distcc %s\" -j %d" physbam-compiler physbam-compile-count))
                                      ((string= physbam-compile-mode "icecream") (format "CXX=\"/opt/icecream/bin/g++\" -j %d" physbam-compile-count))
                                      (physbam-compiler (format "CXX=%s -j %d" physbam-compiler physbam-compile-count))
                                      (t (format "-j %d" physbam-compile-count)))))
      (setq compile-command (format "make -k TYPE=%s %s"
                                physbam-project-type
                                (cond ((string= physbam-compile-mode "distcc") (format "PHYSBAM_CC=\"distcc %s\" -j %d" physbam-compiler physbam-compile-count))
                                      ((string= physbam-compile-mode "icecream") (format "PHYSBAM_CC=\"/opt/icecream/bin/g++\" -j %d" physbam-compile-count))
                                      (t (format "PHYSBAM_CC=%s -j %d" physbam-compiler physbam-compile-count))))))
  (message (format "New compile command is: %s" compile-command))
  (define-key global-map [menu-bar physbam project-dir] '(menu-item (concat  "Project: " physbam-project-directory) physbam-choose-project-directory))
  (define-key global-map [menu-bar physbam output-dir] '(menu-item (concat  "Output: " physbam-output-directory) physbam-choose-output-directory))
  (setq compilation-read-command nil)
  (if (eq write_settings t) (physbam-write-project-settings)))

; commands to store variables
(defun physbam-set-variables-string (variables)
  (physbam-reduce 'concat
          (mapcar (lambda (x) (with-output-to-string (print x) (princ "\n")))
                  (mapcar (lambda (x) `(setq ,x ,(eval x))) variables))))

(defun physbam-write-project-settings ()
  (if (not (string= (replace-regexp-in-string "PhysBAM" "" physbam-project-directory) physbam-project-directory))
      (with-temp-file (format "%s/.emacs_project_config" physbam-project-directory)
        (insert (physbam-set-variables-string '(physbam-project-type
                                                physbam-compile-mode
                                                physbam-compile-count
                                                physbam-output-directory))))))
    
(defun physbam-read-project-settings ()
  (let ((filename (format "%s/.emacs_project_config" physbam-project-directory))) 
    (if (file-exists-p filename)
        (load-file filename))
    (physbam-update-output-directory-menu)))

;#####################################################################
; Menu Bar
;#####################################################################
; Generic helpful scripts
(define-key global-map [menu-bar physbam] (cons "PhysBAM" (make-sparse-keymap "PhysBAM")))
(define-key global-map [menu-bar physbam edit-sconstruct-settings] '("Edit SConstruct" . (lambda () (interactive) (find-file (concat (getenv "PHYSBAM") "/Scripts/scons/SConstruct")))))
(define-key global-map [menu-bar physbam edit-physbam-settings] '("Edit PhysBAM Emacs" . (lambda () (interactive) (find-file (concat (getenv "PHYSBAM") "/Scripts/emacs/physbam.el")))))
(define-key global-map [menu-bar physbam edit-settings] '("Edit .emacs" . (lambda () (interactive) (find-file (concat (getenv "PHYSBAM") "~/.emacs")))))
(define-key global-map [menu-bar physbam sep0] '(menu-item "--single-line"))
(define-key global-map [menu-bar physbam fix-copyright] '("Fix Copyright" . physbam-fix-copyright))
(define-key global-map [menu-bar physbam fix-function-comment] '("Fix Function Comment" . physbam-fix-function-comment))
(define-key global-map [menu-bar physbam update-tags] '("Update Tags" . physbam-update-tags))
(define-key global-map [menu-bar physbam sep1] '(menu-item "--single-line"))
(define-key global-map [menu-bar physbam teeoutput] '(menu-item "Tee Output On Run" (lambda () (interactive) (setq physbam-tee-output (if physbam-tee-output nil t))) :button (:toggle . physbam-tee-output)))
(define-key global-map [menu-bar physbam run] '("Run Program" . physbam-run))
(define-key global-map [menu-bar physbam runviewer] '("Run Viewer" . physbam-run-viewer))
(define-key global-map [menu-bar physbam run-debug] '("Run Debugger" . physbam-run-debug))
(define-key global-map [menu-bar physbam sep1a] '(menu-item "--single-line"))
; Settings menu
(define-key global-map [menu-bar physbam compilesettings] (cons "Compile Settings" (make-sparse-keymap "Compile Settings")))
(define-key global-map [menu-bar physbam compilesettings  mode-distcc] '(menu-item "distcc" (lambda () (interactive) (physbam-set-compile-mode "distcc")) :button (:toggle . (string= physbam-compile-mode "distcc"))))
(define-key global-map [menu-bar physbam compilesettings  mode-icecc] '(menu-item "icecream" (lambda () (interactive) (physbam-set-compile-mode "icecream")) :button (:toggle . (string= physbam-compile-mode "icecream"))))
(define-key global-map [menu-bar physbam compilesettings  mode-single] '(menu-item "single" (lambda () (interactive) (physbam-set-compile-mode "single")) :button (:toggle . (string= physbam-compile-mode "single"))))
(define-key global-map [menu-bar physbam compilesettings sep5] '(menu-item "--single-line"))
(define-key global-map [menu-bar physbam compilesettings j1] '(menu-item "-j 1" (lambda () (interactive) (physbam-set-compile-count 1)) :button (:toggle . (= physbam-compile-count 1))))
(define-key global-map [menu-bar physbam compilesettings j2] '(menu-item "-j 2" (lambda () (interactive) (physbam-set-compile-count 2)) :button (:toggle . (= physbam-compile-count 2))))
(define-key global-map [menu-bar physbam compilesettings j3] '(menu-item "-j 3" (lambda () (interactive) (physbam-set-compile-count 3)) :button (:toggle . (= physbam-compile-count 3))))
(define-key global-map [menu-bar physbam compilesettings j4] '(menu-item "-j 4" (lambda () (interactive) (physbam-set-compile-count 4)) :button (:toggle . (= physbam-compile-count 4))))
(define-key global-map [menu-bar physbam compilesettings j8] '(menu-item "-j 8" (lambda () (interactive) (physbam-set-compile-count 8)) :button (:toggle . (= physbam-compile-count 8))))
(define-key global-map [menu-bar physbam compilesettings j16] '(menu-item "-j 16" (lambda () (interactive) (physbam-set-compile-count 16)) :button (:toggle . (= physbam-compile-count 16))))
(define-key global-map [menu-bar physbam compilesettings j32] '(menu-item "-j 32" (lambda () (interactive) (physbam-set-compile-count 32)) :button (:toggle . (= physbam-compile-count 32))))
(define-key global-map [menu-bar physbam compilesettings j64] '(menu-item "-j 64" (lambda () (interactive) (physbam-set-compile-count 64)) :button (:toggle . (= physbam-compile-count 64))))
(define-key global-map [menu-bar physbam compilesettings sep6] '(menu-item "--single-line"))
(define-key global-map [menu-bar physbam compilesettings g++] '(menu-item "g++" (lambda () (interactive) (physbam-set-compiler "g++")) :button (:toggle . (string= physbam-compiler-id "g++"))))
(define-key global-map [menu-bar physbam compilesettings g++-3.4] '(menu-item "g++ 3.4" (lambda () (interactive) (physbam-set-compiler "gcc-3.4")) :button (:toggle . (string= physbam-compiler-id "gcc-3.4"))))
(define-key global-map [menu-bar physbam compilesettings g++-3.4-64] '(menu-item "g++ 3.4 64" (lambda () (interactive) (physbam-set-compiler "gcc-3.4-64")) :button (:toggle . (string= physbam-compiler-id "gcc-3.4-64"))))
(define-key global-map [menu-bar physbam compilesettings g++-4.0.1] '(menu-item "g++ 3.4" (lambda () (interactive) (physbam-set-compiler "gcc-3.4")) :button (:toggle . (string= physbam-compiler-id "gcc-3.4"))))
(define-key global-map [menu-bar physbam compilesettings g++-4.0.1] '(menu-item "g++ 4.0.1" (lambda () (interactive) (physbam-set-compiler "gcc-4.0.1")) :button (:toggle . (string= physbam-compiler-id "gcc-4.0.1"))))
(define-key global-map [menu-bar physbam compilesettings g++-4.0.1] '(menu-item "g++ 4.1.1" (lambda () (interactive) (physbam-set-compiler "gcc-4.1.1")) :button (:toggle . (string= physbam-compiler-id "gcc-4.1.1"))))
(define-key global-map [menu-bar physbam compilesettings g++-4.0.1-64] '(menu-item "g++ 4.0.1 64" (lambda () (interactive) (physbam-set-compiler "gcc-4.0.1-64")) :button (:toggle . (string= physbam-compiler-id "gcc-4.0.1-64"))))
(define-key global-map [menu-bar physbam compilesettings g++-4.1.1-64] '(menu-item "g++ 4.1.1 64" (lambda () (interactive) (physbam-set-compiler "gcc-4.1.1-64")) :button (:toggle . (string= physbam-compiler-id "gcc-4.1.1-64"))))
(define-key global-map [menu-bar physbam compilesettings g++-3.3.2] '(menu-item "g++ 3.3.2" (lambda () (interactive) (physbam-set-compiler "gcc-3.3.2")) :button (:toggle . (string= physbam-compiler-id "gcc-3.3.2"))))
(define-key global-map [menu-bar physbam compilesettings icc] '(menu-item "icc" (lambda () (interactive) (physbam-set-compiler "icc")) :button (:toggle . (string= physbam-compiler-id "icc"))))
(define-key global-map [menu-bar physbam compilesettings sep7] '(menu-item "--single-line"))
; other compile stuff
(define-key global-map [menu-bar physbam compile] '("Compile Project" . physbam-compile))
(define-key global-map [menu-bar physbam compile-file] '("Compile Current Buffer" . physbam-compile-current-file))
(define-key global-map [menu-bar physbam compile-file-msvc] '("Compile Current Buffer with MSVC" . physbam-compile-current-file-msvc))
(define-key global-map [menu-bar physbam sep1b] '(menu-item "--single-line"))
; Debug or release
(define-key global-map [menu-bar physbam release] '(menu-item "Release" (lambda () (interactive) (physbam-set-project-type "release")) :button (:toggle . (string= physbam-project-type "release"))))
(define-key global-map [menu-bar physbam optdebug] '(menu-item "Optimized Debug" (lambda () (interactive) (physbam-set-project-type "optdebug")) :button (:toggle . (string= physbam-project-type "optdebug"))))
(define-key global-map [menu-bar physbam debug] '(menu-item "Debug" (lambda () (interactive) (physbam-set-project-type "debug")) :button (:toggle . (string= physbam-project-type "debug"))))
(define-key global-map [menu-bar physbam profile] '(menu-item "Profile" (lambda () (interactive) (physbam-set-project-type "profile")) :button (:toggle . (string= physbam-project-type "profile"))))
(define-key global-map [menu-bar physbam sep2] '(menu-item "--single-line"))
(define-key global-map [menu-bar physbam autosetproject] '(menu-item "Set Project From Current" (lambda () (interactive) (physbam-find-project))))

;#####################################################################
; Development related variables
;#####################################################################

(defun physbam-find-project-helper (filename orig-filename depth)
  (cond
   ((file-exists-p (concat filename "SConscript")) (expand-file-name filename))
   ((> depth 10) orig-filename)
   (t (physbam-find-project-helper (concat filename "../") orig-filename (+ depth 1)))))

(defun physbam-find-project ()
  (interactive)
  (setq physbam-project-directory (physbam-find-project-helper default-directory default-directory 0)))

(physbam-find-project)

;(setq physbam-compiler (if (string= (getenv "PLATFORM") "nocona")  "/usr/local/compilers/gcc-3.4-64/bin/g++" "icc"))
(setq physbam-project-type "release")
(setq physbam-compile-count 64)
(setq physbam-tee-output t) ; dump output to buffer on runs
(setq compile-command (format "make -k" physbam-project-directory))
(setq tags-file-name (format "%s/TAGS" (getenv "PHYSBAM")))
(setq physbam-compile-mode "icecream")
; Setup default project to be currenct directory
(setq physbam-output-directory (format "%s/Projects/articulated_rigid_bodies/Blocks/output" (getenv "PHYSBAM")))
; Read project settings
;(physbam-read-project-settings)
; NOTE ALL FUNCTIONS THAT MODIFY STATUS AND SAVE SHOULD BE BELOW ABOVE READ PROJECT SETTINGS
(physbam-set-compiler (if (or (string= (getenv "PLATFORM") "opteron") (string= (getenv "PLATFORM") "nocona"))  "gcc-4.0.1-64" "gcc-4.0.1"))
(physbam-setup-compile-command nil)
(setq truncate-partial-width-windows nil)
(setq compilation-scroll-output t) ; scroll to end by default

;#####################################################################
; Example binds
;#####################################################################
;(global-set-key "\^Xh" 'physbam-header-flip)
;(global-set-key (kbd "<f2>") 'physbam-fix-copyright)
;(global-set-key (kbd "<f3>") 'physbam-fix-function-comment)
;(global-set-key (kbd "<f4>") 'next-error)
;(global-set-key (kbd "<f5>") 'physbam-run)
;(global-set-key (kbd "C-<f5>") 'physbam-run-debug)
;(global-set-key (kbd "<f6>") 'physbam-run-viewer)
;(global-set-key (kbd "<f7>") 'compile)
