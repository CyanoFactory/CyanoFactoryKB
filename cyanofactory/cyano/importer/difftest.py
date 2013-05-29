import diff_match_patch as diff

text1 = """Lorem ipsum dolor sit amet, consectetuer adipiscing elit. Integer
eu lacus accumsan arcu fermentum euismod. Donec pulvinar porttitor
tellus. Aliquam venenatis. Donec facilisis pharetra tortor.  In nec
mauris eget magna consequat convallis. Nam sed sem vitae odio
pellentesque interdum. Sed consequat viverra nisl. Suspendisse arcu
metus, blandit quis, rhoncus ac, pharetra eget, velit. Mauris
urna. Morbi nonummy molestie orci. Praesent nisi elit, fringilla ac,
suscipit non, tristique vel, mauris. Curabitur vel lorem id nisl porta
adipiscing. Suspendisse eu lectus. In nunc. Duis vulputate tristique
enim. Donec quis lectus a justo imperdiet tempus."""
text1_lines = text1.splitlines()

text2 = """Lorem ipsum dolor sit amet, consectetuer adipiscing elit. Integer
eu lacus accumsan arcu fermentum euismod. Donec pulvinar, porttitor
tellus. Aliquam venenatis. Donec facilisis pharetra tortor. In nec
mauris eget magna consequat convallis. Nam cras vitae mi vitae odio
pellentesque interdum. Sed consequat viverra nisl. Suspendisse arcu
metus, blandit quis, rhoncus ac, pharetra eget, velit. Mauris
urna. Morbi nonummy molestie orci. Praesent nisi elit, fringilla ac,
suscipit non, tristique vel, mauris. Curabitur vel lorem id nisl porta
adipiscing. Duis vulputate tristique enim. Donec quis lectus a justo
imperdiet tempus. Suspendisse eu lectus. In nunc. """
text2_lines = text2.splitlines()

d = diff.diff_match_patch()
di = d.diff_main(text1.replace("\n", " "), text2.replace("\n", " "))
print di

f = open("diff.html", "w")
f.write("<table>")
f.write("<tr>")
f.write("<td>")
for i, t in di:
    if i == -1:
        f.write('<span style="background: red">' + t + '</span>')
    elif i == 0:
        f.write(t)
    else:
        pass
f.write("</td>")
f.write("<td>")
for i, t in di:
    if i == -1:
        pass
    elif i == 0:
        f.write(t)
    else:
        f.write('<span style="background: green">' + t + '</span>')
f.write("</td>")
f.write("</tr>")
f.write("</table>")
f.close()