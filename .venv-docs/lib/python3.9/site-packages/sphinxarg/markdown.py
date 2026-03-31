try:
    from commonmark import Parser
except ImportError:
    from CommonMark import Parser  # >= 0.5.6
try:
    from commonmark.node import Node
except ImportError:
    from CommonMark.node import Node
from docutils import nodes
from docutils.utils.code_analyzer import Lexer


def custom_walker(node, space=''):
    """
    A convenience function to ease debugging. It will print the node structure that's returned from CommonMark

    The usage would be something like:

    >>> content = Parser().parse('Some big text block\n===================\n\nwith content\n')
    >>> custom_walker(content)
    document
        heading
            text    Some big text block
        paragraph
            text    with content

    Spaces are used to convey nesting
    """
    txt = ''
    try:
        txt = node.literal
    except Exception:
        pass

    if txt is None or txt == '':
        print(f'{space}{node.t}')
    else:
        print(f'{space}{node.t}\t{txt}')

    cur = node.first_child
    if cur:
        while cur is not None:
            custom_walker(cur, space + '    ')
            cur = cur.nxt


def paragraph(node):
    """
    Process a paragraph, which includes all content under it
    """
    text = ''
    if node.string_content is not None:
        text = node.string_content
    o = nodes.paragraph('', ' '.join(text))
    o.line = node.sourcepos[0][0]
    for n in markdown(node):
        o.append(n)

    return o


def text(node):
    """
    Text in a paragraph
    """
    return nodes.Text(node.literal)


def hardbreak(node):
    """
    A <br /> in html or "\n" in ascii
    """
    return nodes.Text('\n')


def softbreak(node):
    """
    A line ending or space.
    """
    return nodes.Text('\n')


def reference(node):
    """
    A hyperlink. Note that alt text doesn't work, since there's no apparent way to do that in docutils
    """
    o = nodes.reference()
    o['refuri'] = node.destination
    if node.title:
        o['name'] = node.title
    for n in markdown(node):
        o += n
    return o


def emphasis(node):
    """
    An italicized section
    """
    o = nodes.emphasis()
    for n in markdown(node):
        o += n
    return o


def strong(node):
    """
    A bolded section
    """
    o = nodes.strong()
    for n in markdown(node):
        o += n
    return o


def literal(node):
    """
    Inline code
    """
    rendered = []
    try:
        if node.info is not None:
            l = Lexer(node.literal, node.info, tokennames="long")
            for _ in l:
                rendered.append(node.inline(classes=_[0], text=_[1]))
    except Exception:
        pass

    classes = ['code']
    if node.info is not None:
        classes.append(node.info)
    if len(rendered) > 0:
        o = nodes.literal(classes=classes)
        for element in rendered:
            o += element
    else:
        o = nodes.literal(text=node.literal, classes=classes)

    for n in markdown(node):
        o += n
    return o


def literal_block(node):
    """
    A block of code
    """
    rendered = []
    try:
        if node.info is not None:
            l = Lexer(node.literal, node.info, tokennames="long")
            for _ in l:
                rendered.append(node.inline(classes=_[0], text=_[1]))
    except Exception:
        pass

    classes = ['code']
    if node.info is not None:
        classes.append(node.info)
    if len(rendered) > 0:
        o = nodes.literal_block(classes=classes)
        for element in rendered:
            o += element
    else:
        o = nodes.literal_block(text=node.literal, classes=classes)

    o.line = node.sourcepos[0][0]
    for n in markdown(node):
        o += n
    return o


def raw(node):
    """
    Add some raw html (possibly as a block)
    """
    o = nodes.raw(node.literal, node.literal, format='html')
    if node.sourcepos is not None:
        o.line = node.sourcepos[0][0]
    for n in markdown(node):
        o += n
    return o


def transition(node):
    """
    An <hr> tag in html. This has no children
    """
    return nodes.transition()


def title(node):
    """
    A title node. It has no children
    """
    return nodes.title(node.first_child.literal, node.first_child.literal)


def section(node):
    """
    A section in reStructuredText, which needs a title (the first child)
    This is a custom type
    """
    title = ''  # All sections need an id
    if node.first_child is not None:
        if node.first_child.t == 'heading':
            title = node.first_child.first_child.literal
    o = nodes.section(ids=[title], names=[title])
    for n in markdown(node):
        o += n
    return o


def block_quote(node):
    """
    A block quote
    """
    o = nodes.block_quote()
    o.line = node.sourcepos[0][0]
    for n in markdown(node):
        o += n
    return o


def image(node):
    """
    An image element

    The first child is the alt text. reStructuredText can't handle titles
    """
    o = nodes.image(uri=node.destination)
    if node.first_child is not None:
        o['alt'] = node.first_child.literal
    return o


def list_item(node):
    """
    An item in a list
    """
    o = nodes.list_item()
    for n in markdown(node):
        o += n
    return o


def list_node(node):
    """
    A list (numbered or not)
    For numbered lists, the suffix is only rendered as . in html
    """
    if node.list_data['type'] == 'bullet':
        o = nodes.bullet_list(bullet=node.list_data['bullet_char'])
    else:
        o = nodes.enumerated_list(
            suffix=node.list_data['delimiter'],
            enumtype='arabic',
            start=node.list_data['start'],
        )
    for n in markdown(node):
        o += n
    return o


def markdown(node):
    """
    Returns a list of nodes, containing CommonMark nodes converted to docutils nodes
    """
    cur = node.first_child

    # Go into each child, in turn
    output = []
    while cur is not None:
        t = cur.t
        if t == 'paragraph':
            output.append(paragraph(cur))
        elif t == 'text':
            output.append(text(cur))
        elif t == 'softbreak':
            output.append(softbreak(cur))
        elif t == 'linebreak':
            output.append(hardbreak(cur))
        elif t == 'link':
            output.append(reference(cur))
        elif t == 'heading':
            output.append(title(cur))
        elif t == 'emph':
            output.append(emphasis(cur))
        elif t == 'strong':
            output.append(strong(cur))
        elif t == 'code':
            output.append(literal(cur))
        elif t == 'code_block':
            output.append(literal_block(cur))
        elif t == 'html_inline' or t == 'html_block':
            output.append(raw(cur))
        elif t == 'block_quote':
            output.append(block_quote(cur))
        elif t == 'thematic_break':
            output.append(transition(cur))
        elif t == 'image':
            output.append(image(cur))
        elif t == 'list':
            output.append(list_node(cur))
        elif t == 'item':
            output.append(list_item(cur))
        elif t == 'MDsection':
            output.append(section(cur))
        else:
            print(f'Received unhandled type: {t}. Full print of node:')
            cur.pretty()

        cur = cur.nxt

    return output


def finalize_section(section):
    """
    Correct the nxt and parent for each child
    """
    cur = section.first_child
    last = section.last_child
    if last is not None:
        last.nxt = None

    while cur is not None:
        cur.parent = section
        cur = cur.nxt


def nest_sections(block, level=1):
    """
    Sections aren't handled by CommonMark at the moment.
    This function adds sections to a block of nodes.
    'title' nodes with an assigned level below 'level' will be put in a child section.
    If there are no child nodes with titles of level 'level' then nothing is done
    """
    cur = block.first_child
    if cur is not None:
        children = []
        # Do we need to do anything?
        nest = False
        while cur is not None:
            if cur.t == 'heading' and cur.level == level:
                nest = True
                break
            cur = cur.nxt
        if not nest:
            return

        section = Node('MDsection', 0)
        section.parent = block
        cur = block.first_child
        while cur is not None:
            if cur.t == 'heading' and cur.level == level:
                # Found a split point, flush the last section if needed
                if section.first_child is not None:
                    finalize_section(section)
                    children.append(section)
                    section = Node('MDsection', 0)
            nxt = cur.nxt
            # Avoid adding sections without titles at the start
            if section.first_child is None:
                if cur.t == 'heading' and cur.level == level:
                    section.append_child(cur)
                else:
                    children.append(cur)
            else:
                section.append_child(cur)
            cur = nxt

        # If there's only 1 child then don't bother
        if section.first_child is not None:
            finalize_section(section)
            children.append(section)

        block.first_child = None
        block.last_child = None
        next_level = level + 1
        for child in children:
            # Handle nesting
            if child.t == 'MDsection':
                nest_sections(child, level=next_level)

            # Append
            if block.first_child is None:
                block.first_child = child
            else:
                block.last_child.nxt = child
            child.parent = block
            child.nxt = None
            child.prev = block.last_child
            block.last_child = child


def parse_markdown_block(text):
    """
    Parses a block of text, returning a list of docutils nodes

    >>> parse_markdown_block("Some\n====\n\nblock of text\n\nHeader\n======\n\nblah\n")
    []
    """
    block = Parser().parse(text)
    # CommonMark can't nest sections, so do it manually
    nest_sections(block)

    return markdown(block)
