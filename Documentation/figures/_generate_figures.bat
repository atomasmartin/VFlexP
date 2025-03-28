for %%F in (*.tex) do (
  if exist "%%~dpnF.pdf" (
    for %%P in ("%%~dpnF.pdf") do (
      if "%%~tF" gtr "%%~tP" (
        pdflatex "%%F"
      )
    )
  ) else (
    pdflatex "%%F"
  )
)

FOR /R %%F IN (*.pdf) DO (
  IF EXIST "%%~dpnF.png" (
    FOR /F %%D IN ('DIR /TC "%%~dpnF.png" ^| FIND /I "%%~nF.png"') DO (
      IF "%%~tF" GTR "%%D" (
        magick convert -density 300 "%%F" "%%~dpnF.png"
      )
    )
  ) ELSE (
    magick convert -density 300 "%%F" "%%~dpnF.png"
  )
)


