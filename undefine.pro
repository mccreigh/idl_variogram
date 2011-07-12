;+
; NAME:
;       UNDEFINE
;
; PURPOSE:
;       The purpose of this program is to delete or undefine
;       an IDL program variable from within an IDL program or
;       at the IDL command line. It is a more powerful DELVAR.
;
; AUTHOR:
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       1642 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com
;
; CATEGORY:
;       Utilities.
;
; CALLING SEQUENCE:
;       UNDEFINE, variable
;
; REQUIRED INPUTS:
;       variable: The variable to be deleted. Up to 10 variables may be specified as arguments.
;
; SIDE EFFECTS:
;       The variable no longer exists.
;
; EXAMPLE:
;       To delete the variable "info", type:
;
;        IDL> Undefine, info
;
; MODIFICATION HISTORY:
;       Written by David Fanning, 8 June 97, from an original program
;       given to me by Andrew Cool, DSTO, Adelaide, Australia.
;       Simplified program so you can pass it an undefined variable. :-) 17 May 2000. DWF
;       Simplified it even more by removing the unnecessary SIZE function. 28 June 2002. DWF.
;       Added capability to delete up to 10 variables at suggestion of Craig Markwardt. 10 Jan 2008. DWF.
;-
PRO UNDEFINE, var0, var1, var2, var3, var4, var5, var6, var7, var8, var9
   var0 = 0 & dummy = Temporary(var0)
   var1 = 0 & dummy = Temporary(var1)
   var2 = 0 & dummy = Temporary(var2)
   var3 = 0 & dummy = Temporary(var3)
   var4 = 0 & dummy = Temporary(var4)
   var5 = 0 & dummy = Temporary(var5)
   var6 = 0 & dummy = Temporary(var6)
   var7 = 0 & dummy = Temporary(var7)
   var8 = 0 & dummy = Temporary(var8)
   var9 = 0 & dummy = Temporary(var9)
END