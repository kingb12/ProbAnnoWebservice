#! /usr/bin/env python
import cgi
import cgitb
cgitb.enable()
#print """Content-type: text/html
print """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<title>ProbAnno</title>
</head>
<body>
<form method="post" action="ProbAnno.cgi" enctype="multipart/form-data">

  <h3> ProbAnno </h3>
<p> This is a test!</p>
<p>Your results are <a href="http://www.google.com">here</a>
    <p>Uniprot Proteome ID (e.g. UP000018540): <input type="text" name="UniprotID"/>
<br><i>
Browse available proteomes at <a target="_blank" href="http://www.uniprot.org/proteomes/">UniProt</a></i>
<p>OR upload a proteome fasta file: <input type="file" name="uploaded_file"  size="30" maxlength="200" />
</body>
</html>
""" 
